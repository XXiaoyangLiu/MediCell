#!/usr/bin/env python3
"""MediCell pipeline driver: chains the M01-M08 scripts end to end.

All paths are passed to the pipeline scripts via MC_* environment variables.
Intermediate and final results are written under --output-dir with fixed
step-based names. Every step runs sequentially and fails fast.
Use --cell-type-twas (the default) to run the cell-type TWAS pipeline, or --steps STEP [STEP ...]
for a subset in pipeline order (space or comma separated).

Canonical step order:
    M01 M02exp M02cov M03 M04db M05align M05metal M05filter M06ld M06cov M07 M08
Old step names still work as aliases: M06 (database), M04align/M04metal/M04filter
(GWAS), M05ld/M05cov (LD). Historical output directory names are unchanged.

GWAS steps (M05align/M05metal/M05filter) require an explicit mode:
  --gwas-profile   bit-exact reproduction of the historical SLE implementation,
                   including its known statistical errors (legacy only)
  --gwas-manifest  corrected standard path from an explicit study manifest
The active mode is recorded in M04_gwas/gwas_run.json so M05metal/M05filter/M07
can resume it when re-run in the same output directory.

Usage example:
    python3 medicell.py --output-dir ./run_out --cell-type-twas \
        --cell-types B,CD4T,CD8T,NK,Mono,DC --peer-nk 2 --peer-keep 2 \
        --training-population EUR --twas-model-set unfiltered \
        --bulk-tpm bulk_TPM.txt --enigma-ref sc_ref.matrx \
        --cibersortx-frax CIBERSORTx_frax.txt --covariates covs.txt \
        --snp-annot-pattern "snp_annot/snp_annot.chr{chr}.txt" \
        --genotype-pattern "genotype/genotype.chr{chr}.txt" \
        --gene-annot gene_annot.parsed.txt --snp-annot-file snp_annotation.txt \
        --gwas-dir gwas/ --gwas-profile sle-legacy-pbmc \
        --g1000-vcf-pattern "ALL.chr{chr}.phase3.vcf.gz" \
        --g1000-samples EUR_EAS_samples.txt --gene-pos gene_pos_grch37.82.csv
"""
import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import tempfile
import threading

REPO = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(REPO, "script")

ALL_STEPS = [
    "inferExp", "normExp", "peerCov", "trainModel", "buildDB",
    "alignGWAS", "metaGWAS", "filterGWAS",
    "extractLD", "buildLD", "runTWAS", "summarizeTWAS",
]

STEP_DESCRIPTIONS = {
    "inferExp": "deconvolve bulk TPM into cell-type expression (ENIGMA)",
    "normExp": "filter expressed genes and INT-transform per cell type",
    "peerCov": "estimate PEER hidden factors and format covariates",
    "trainModel": "train cell-type elastic-net GReX models",
    "buildDB": "build PrediXcan SQLite weight databases",
    "alignGWAS": "align GWAS summary statistics to the SNP panel",
    "metaGWAS": "run METAL GWAS meta-analysis",
    "filterGWAS": "filter meta-analyzed GWAS variants",
    "extractLD": "extract and align 1000G reference dosages",
    "buildLD": "compute cell-type LD covariance matrices",
    "runTWAS": "run S-PrediXcan cell-type TWAS",
    "summarizeTWAS": "summarize TWAS results with FDR and Bonferroni",
}

# Historical M0* step names keep working as aliases.
STEP_ALIASES = {
    "M01": "inferExp",
    "M02exp": "normExp",
    "M02cov": "peerCov",
    "M03": "trainModel",
    "M04db": "buildDB",
    "M06": "buildDB",
    "M05align": "alignGWAS",
    "M04align": "alignGWAS",
    "M05metal": "metaGWAS",
    "M04metal": "metaGWAS",
    "M05filter": "filterGWAS",
    "M04filter": "filterGWAS",
    "M06ld": "extractLD",
    "M05ld": "extractLD",
    "M06cov": "buildLD",
    "M05cov": "buildLD",
    "M07": "runTWAS",
    "M08": "summarizeTWAS",
}

SCRIPT_FILES = {
    "inferExp": "M01_ENIGMA.R",
    "normExp": "M02_DataProcess_exp_fordb.R",
    "peerCov": "M02_DataProcess_covariance_for_db.R",
    "trainModel": "M03_ModelTraining.R",
    "buildDB": "M04_Makedb.R",
    "alignGWAS": "M05_DataProcess_GWAS_align4METAL.R",
    "metaGWAS": "M05_DataProcess_GWAS_metal.sh",
    "filterGWAS": "M05_DataProcess_GWAS_metal_filter.R",
    "extractLD": "M06_LD.sh",
    "buildLD": "M06_Makeld.sh",
    "runTWAS": "M07_TWAS.sh",
    "summarizeTWAS": "M08_ResultStatistic.R",
}

LINCS_STEPS = ["prepLincs", "enrichLincs", "scoreLincs"]
LINCS_STEP_ALIASES = {
    "L01": "prepLincs",
    "L02": "enrichLincs",
    "L03": "scoreLincs",
}
LINCS_SCRIPT_FILES = {
    "prepLincs": "L01_lincs_prepare.R",
    "enrichLincs": "L02_lincs_enrichment.R",
    "scoreLincs": "L03_lincs_drug_score.sh",
}

GWAS_STEPS = {"alignGWAS", "metaGWAS", "filterGWAS"}
LEGACY_GWAS_PREFIX = "SLE_metal_096_917_183"
LEGACY_GWAS_PROFILES = ("sle-legacy-pbmc", "sle-legacy-kidney")
SAFE_NAME = r"[A-Za-z0-9][A-Za-z0-9_.\-]*"


def write_metadata(path, data):
    with tempfile.NamedTemporaryFile(mode="w", dir=os.path.dirname(path),
                                     prefix=".gwas-run-", delete=False) as f:
        temporary = f.name
        try:
            json.dump(data, f, indent=2, sort_keys=True)
            f.flush()
            os.fsync(f.fileno())
        except BaseException:
            os.unlink(temporary)
            raise
    try:
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read_metadata(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="MediCell pipeline driver (use `medicell.py lincs --help` for LINCS therapeutic prioritization)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--alpha", type=float, default=0.1,
                   help="ENIGMA alpha (default 0.1)")
    p.add_argument("--tpm-min", type=float, default=0.1,
                   help="normExp gene filter: min inferred TPM (default 0.1)")
    p.add_argument("--min-pct", type=float, default=0.01,
                   help="normExp gene filter: min sample fraction, 0-1 (default 0.01)")
    p.add_argument("--min-samples", type=int, default=5,
                   help="normExp gene filter: min sample count floor (default 5)")
    p.add_argument("--cat-covs", default="",
                   help="peerCov: comma-separated categorical covariate columns to dummy-encode")
    p.add_argument("--min-level-count", type=int, default=5,
                   help="peerCov: drop categorical levels with fewer samples (default 5)")
    p.add_argument("--peer-nk", type=int, default=None,
                   help="PEER hidden factors; required explicitly for peerCov")
    p.add_argument("--peer-keep", type=int, default=None,
                   help="PEER factors kept; required explicitly for peerCov")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--cell-type-twas", action="store_true",
                      help="run the cell-type TWAS pipeline (default when no step selector is given)")
    mode.add_argument("--all", dest="cell_type_twas", action="store_true",
                      help=argparse.SUPPRESS)
    mode.add_argument("--steps", nargs="+", metavar="STEP",
                      help="run only these steps, separated by spaces or commas, "
                           "in pipeline order: " + ",".join(ALL_STEPS)
                           + " (old M01-M08 names remain valid aliases)")
    p.add_argument("--chromosomes", default=",".join(str(c) for c in range(1, 23)),
                   help="comma-separated chromosome list (default 1-22)")
    p.add_argument("--cell-types", required=True,
                   help="space/comma-separated cell types matching the reference and fractions")

    data = p.add_argument_group("data inputs")
    data.add_argument("--bulk-tpm", default="",
                      help="bulk TPM expression matrix (genes x samples, tab-separated)")
    data.add_argument("--enigma-ref", default="",
                      help="single-cell reference matrix (TPM, aggregated by cell type)")
    data.add_argument("--cibersortx-frax", default="",
                      help="CIBERSORTx cell-fraction output estimated externally and passed to ENIGMA")
    data.add_argument("--covariates", default="",
                      help="covariates txt (required by peerCov auto-detects "
                           "samples-x-covariates and legacy covariates-x-samples layouts; "
                           "see README)")
    data.add_argument("--snp-annot-pattern", default="",
                      help="per-chromosome SNP annotation, {chr} placeholder, e.g. snp_annot.chr{chr}.txt")
    data.add_argument("--genotype-pattern", default="",
                      help="per-chromosome genotype dosage, {chr} placeholder, e.g. genotype.chr{chr}.txt")
    data.add_argument("--gene-annot", default="",
                      help="gene annotation file (parsed, e.g. gene_annot.parsed.txt)")
    data.add_argument("--snp-annot-file", default="",
                      help="single-file genotype panel SNP annotation (used by the GWAS align step)")
    data.add_argument("--gwas-dir", default="",
                      help="directory holding the historical SLE GWAS files "
                           "(used by alignGWAS)")
    data.add_argument("--g1000-vcf-pattern", default="",
                      help="1000 Genomes Phase 3 VCF pattern with {chr} placeholder "
                           "(not needed when --ld-dosage-dir is given)")
    data.add_argument("--g1000-samples", default="",
                      help="1000G sample list matching the GWAS ancestry (e.g., 1000G EUR samples for EUR GWAS)")
    data.add_argument("--gene-pos", default="",
                      help="gene position annotation CSV (GRCh37, gene_id/gene_name columns)")
    data.add_argument("--ld-dosage-dir", default="",
                      help="pre-extracted 1000G dosage directory "
                           "(skips the extractLD extraction block)")
    data.add_argument("--gwas-summary", default="",
                      help="pre-QC'd GWAS summary for runTWAS (MarkerName, Allele1, Allele2, Effect, StdErr, P-value; skips alignGWAS/metaGWAS/filterGWAS with --cell-type-twas)")

    gen = p.add_argument_group("input contracts and generalization")
    gen.add_argument("--sample-id-policy", choices=["exact", "legacy-rnames"], default=None,
                     help="sample ID matching between expression and covariate/genotype "
                          "files (default exact; legacy-rnames must be selected explicitly)")
    gen.add_argument("--cov-layout", choices=["auto", "samples", "covariates"], default="auto",
                     help="peerCov covariates file layout (default auto)")
    gen.add_argument("--gene-annot-format", choices=["standard", "legacy-name-is-id"],
                     default=None,
                     help="gene annotation schema (default standard; legacy-name-is-id must be selected explicitly)")
    gen.add_argument("--genome-build", choices=["GRCh37", "GRCh38"], default=None,
                     help="genome build declared for this cohort's coordinates")
    gen.add_argument("--training-population", default="unspecified",
                     help="optional population label stored with the models (default unspecified)")
    gen.add_argument("--maf", type=float, default=0.01,
                     help="minor allele frequency threshold for trainModel genotype variants in [0, 0.5] (default 0.01)")
    gen.add_argument("--maf-policy", choices=["filter", "verify"], default="filter",
                     help="filter: drop variants with MAF < --maf; verify: fail if input genotype violates --maf (default filter)")
    gen.add_argument("--gwas-population", default=None,
                     help="GWAS ancestry label (e.g. EUR, EAS, EUR_EAS); defaults --ld-panel-id when omitted")
    gen.add_argument("--ld-panel-id", default=None,
                     help="LD panel id matched to GWAS ancestry (defaults to --gwas-population, else EUR_EAS)")
    gen.add_argument("--twas-model-set", choices=["unfiltered", "filtered-signif"], default="filtered-signif",
                     help="weight database for runTWAS (default filtered-signif: R2 > 0.01 and CV p < 0.05; buildDB always keeps *_models_unfiltered.db for buildLD)")
    gen.add_argument("--cores", type=int, default=None,
                     help="unified parallel cores/workers across all pipeline steps (default 4)")
    gen.add_argument("--seed", type=int, default=None,
                     help="explicit trainModel training seed; forces serial inner training "
                          "(--training-cores 1)")
    gen.add_argument("--training-cores", type=int, default=None,
                     help="inner trainModel cores per (cell_type, chr) job (defaults from --cores; forced to 1 with --seed)")

    gwas = p.add_argument_group("GWAS mode (required when running alignGWAS/metaGWAS/filterGWAS)")
    gwas_mode = gwas.add_mutually_exclusive_group()
    gwas_mode.add_argument("--gwas-profile", choices=list(LEGACY_GWAS_PROFILES), default=None,
                           help="historical SLE reproduction profile; reproduces the "
                                "original implementation including its known statistical "
                                "errors - for reproduction only")
    gwas_mode.add_argument("--gwas-manifest", default=None, metavar="FILE",
                           help="TSV manifest of GWAS studies for the corrected "
                                "standard path (see README for the column contract)")
    gwas.add_argument("--gwas-filter-policy", choices=["none", "legacy"], default=None,
                      help="post-meta filtering (default: legacy for --gwas-profile, "
                           "none for --gwas-manifest)")
    gwas.add_argument("--skip-gwas-meta", action="store_true",
                      help="skip optional metaGWAS and filterGWAS when running alignGWAS on a single pre-QC'd GWAS study in --gwas-manifest")

    env = p.add_argument_group("environment and tools")
    env.add_argument("--env-base", default="",
                     help="portable environment root (unpacked conda-pack bundle); "
                          "when set, env binaries are used directly instead of conda run")
    env.add_argument("--r-env", default="R4.4.3", help="conda env with R 4.4.3 + ENIGMA")
    env.add_argument("--peer-env", default="r3.6.3", help="conda env with R 3.6.3 + PEER")
    env.add_argument("--py-env", default="imlatools", help="conda env with MetaXcan runtime")
    env.add_argument("--metaxcan-dir", default=os.path.join(REPO, "tools", "MetaXcan"))
    env.add_argument("--metal-bin", default=os.path.join(REPO, "tools", "METAL", "metal"))
    env.add_argument("--bcftools", default="bcftools", help="bcftools executable")
    env.add_argument("--bcf-plugins", default="", help="bcftools plugins directory (empty = auto)")
    env.add_argument("--gwas-id", default=None,
                     help="trait label used in standard-mode output names and reports "
                          "(the legacy profiles keep their historical output names)")
    step_flags = p.add_argument_group("step shortcuts (alternative to --cell-type-twas / --steps; e.g. medicell.py --inferExp ...)")
    for step in ALL_STEPS:
        step_flags.add_argument("--" + step, dest="flag_" + step, action="store_true",
                                help="run " + step + ": " + STEP_DESCRIPTIONS[step])

    args = p.parse_args(argv)

    # -- steps: normalize shortcuts/aliases, validate, deduplicate, pipeline order --
    flagged_steps = [step for step in ALL_STEPS if getattr(args, "flag_" + step)]
    if flagged_steps and (args.cell_type_twas or args.steps is not None):
        p.error("step shortcut flags (--" + ",--".join(flagged_steps) + ") cannot be combined with --cell-type-twas or --steps")

    args.step_aliases = []
    if flagged_steps:
        args.steps = flagged_steps
    elif args.steps is None:
        if args.gwas_summary:
            args.steps = [step for step in ALL_STEPS if step not in GWAS_STEPS]
        elif args.skip_gwas_meta:
            args.steps = [step for step in ALL_STEPS if step not in ("metaGWAS", "filterGWAS")]
        else:
            args.steps = ALL_STEPS.copy()
    else:
        selected = []
        for value in args.steps:
            for part in value.split(","):
                if not part.strip():
                    p.error("--steps cannot contain empty steps; choose from: "
                            + ",".join(ALL_STEPS))
                selected.extend(part.split())
        canonical = []
        unknown = []
        for s in selected:
            c = STEP_ALIASES.get(s, s)
            if c not in ALL_STEPS:
                unknown.append(s)
                continue
            if s != c and (s, c) not in args.step_aliases:
                args.step_aliases.append((s, c))
            if c not in canonical:
                canonical.append(c)
        if unknown:
            p.error("unknown steps: " + ",".join(sorted(set(unknown)))
                    + "; choose from: " + ",".join(ALL_STEPS))
        args.steps = [step for step in ALL_STEPS if step in canonical]

    # -- numeric and range validation --
    if not math.isfinite(args.alpha) or args.alpha <= 0:
        p.error("--alpha must be a positive finite number")
    if not (0 <= args.min_pct <= 1):
        p.error("--min-pct is a sample fraction between 0 and 1")
    if not math.isfinite(args.tpm_min) or args.tpm_min < 0:
        p.error("--tpm-min must be a finite non-negative number")
    if args.min_samples < 0:
        p.error("--min-samples must be non-negative")
    if args.min_level_count < 0:
        p.error("--min-level-count must be non-negative")
    if not math.isfinite(args.maf) or not 0 <= args.maf <= 0.5:
        p.error("--maf must be a finite number in [0, 0.5]")
    if args.seed is not None and not 0 <= args.seed <= 2147483647:
        p.error("--seed must be an integer in 0..2147483647")
    if args.cores is not None and args.cores < 1:
        p.error("--cores must be >= 1")
    if args.training_cores is not None and args.training_cores < 1:
        p.error("--training-cores must be >= 1")
    if args.seed is not None and args.training_cores not in (None, 1):
        p.error("--seed requires serial inner training; set --training-cores 1 or omit it")
    args.cores_resolved = (args.cores if args.cores is not None
                           else (args.training_cores if args.training_cores is not None else 4))
    chrs = []
    for token in args.chromosomes.split(","):
        token = token.strip()
        if not token:
            p.error("--chromosomes cannot contain empty values")
        if not re.fullmatch(r"\d+", token) or not 1 <= int(token) <= 22:
            p.error("invalid chromosome %r (must be 1-22)" % token)
        if int(token) not in chrs:
            chrs.append(int(token))
    args.chromosomes = chrs
    for name, value in (("--gwas-id", args.gwas_id),
                        ("--gwas-population", args.gwas_population),
                        ("--ld-panel-id", args.ld_panel_id)):
        if value is not None and not re.fullmatch(SAFE_NAME, value):
            p.error(name + " must be a safe output label ([A-Za-z0-9_.-])")

    # -- explicit analysis configuration --
    cell_types = args.cell_types.replace(",", " ").split()
    if not cell_types:
        p.error("--cell-types cannot be empty")
    if any(not re.fullmatch(SAFE_NAME, ct) for ct in cell_types):
        p.error("cell type names must be safe output labels ([A-Za-z0-9_.-])")
    if len(set(cell_types)) != len(cell_types):
        p.error("--cell-types contains duplicates")
    args.cell_types = cell_types

    peer_nk, peer_keep = args.peer_nk, args.peer_keep
    if "peerCov" in args.steps and (peer_nk is None or peer_keep is None):
        p.error("--peer-nk and --peer-keep are required for peerCov")
    if peer_nk is not None and peer_nk < 0:
        p.error("--peer-nk must be non-negative")
    if peer_keep is not None and (peer_keep < 0 or (peer_nk is not None and peer_keep > peer_nk)):
        p.error("--peer-keep must be between 0 and --peer-nk")
    args.peer_nk_resolved, args.peer_keep_resolved = peer_nk, peer_keep
    args.sample_id_policy = args.sample_id_policy or "exact"
    args.gene_annot_format = args.gene_annot_format or "standard"
    if any(step in args.steps for step in ("buildDB", "runTWAS")) and not args.twas_model_set:
        p.error("--twas-model-set (unfiltered or filtered-signif) is required for M04db and M07")

    # -- per-step required inputs --
    required_inputs = {
        "inferExp": ("bulk_tpm", "enigma_ref", "cibersortx_frax"),
        "peerCov": ("covariates",),
        "trainModel": ("snp_annot_pattern", "genotype_pattern", "gene_annot"),
        "alignGWAS": (("snp_annot_file",) if args.gwas_manifest else
                     ("snp_annot_file", "gwas_dir")),
        "extractLD": (("snp_annot_pattern",) if args.ld_dosage_dir else
                  ("snp_annot_pattern", "g1000_vcf_pattern", "g1000_samples")),
        "summarizeTWAS": ("gene_pos",),
    }
    missing = sorted({"--" + name.replace("_", "-")
                      for step in args.steps for name in required_inputs.get(step, ())
                      if not getattr(args, name)})
    if missing:
        p.error("selected steps require: " + ", ".join(missing))

    # -- GWAS mode: external summary, single-study aligned, multi-study meta, or recorded prior mode --
    if args.gwas_summary:
        if args.gwas_profile or args.gwas_manifest:
            p.error("--gwas-summary cannot be combined with --gwas-profile or --gwas-manifest")
        if args.skip_gwas_meta:
            p.error("--skip-gwas-meta cannot be combined with --gwas-summary")
        if set(args.steps) & GWAS_STEPS:
            p.error("--gwas-summary cannot be combined with alignGWAS, metaGWAS, or filterGWAS")
        if not os.path.isfile(args.gwas_summary) or os.path.getsize(args.gwas_summary) == 0:
            p.error("--gwas-summary file not found or empty: " + args.gwas_summary)
        import gzip as _gzip
        opener = _gzip.open if args.gwas_summary.endswith(".gz") else open
        with opener(args.gwas_summary, "rt") as _gf:
            _hdr = _gf.readline().rstrip("\r\n").split()
        _req_cols = {"MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value"}
        if not (_req_cols <= set(_hdr)):
            p.error("--gwas-summary header must contain columns: " + ", ".join(sorted(_req_cols)))
    if args.skip_gwas_meta:
        if args.gwas_profile:
            p.error("--skip-gwas-meta applies only to single-study --gwas-manifest (not --gwas-profile)")
        if set(args.steps) & {"metaGWAS", "filterGWAS"}:
            p.error("--skip-gwas-meta cannot be combined with metaGWAS or filterGWAS")
    gwas_dir = os.path.join(os.path.abspath(args.output_dir), "M04_gwas")
    run_meta = os.path.join(gwas_dir, "gwas_run.json")
    profile, manifest = args.gwas_profile, args.gwas_manifest
    needs_mode = bool(set(args.steps) & GWAS_STEPS) or \
        ("runTWAS" in args.steps and not args.gwas_summary)
    needs_summary_id = "summarizeTWAS" in args.steps
    args.gwas_mode_recovered = False
    saved = None
    requested_profile, requested_manifest = profile, manifest
    if needs_mode or (needs_summary_id and not args.gwas_summary):
        saved = read_metadata(run_meta)
        if saved is not None:
            if not isinstance(saved, dict):
                p.error("Invalid GWAS mode record: " + run_meta)
            saved_id = saved.get("gwas_id")
            if not isinstance(saved_id, str) or not re.fullmatch(SAFE_NAME, saved_id):
                p.error("Invalid GWAS identity record: " + run_meta)
            if args.gwas_id is not None and args.gwas_id != saved_id:
                p.error("--gwas-id conflicts with recorded GWAS inputs; use a new output directory")
            if needs_mode:
                profile, manifest = saved.get("profile"), saved.get("manifest")
                if (bool(profile) == bool(manifest)
                        or (profile and profile not in LEGACY_GWAS_PROFILES)
                        or (manifest and (not isinstance(manifest, str)
                                          or not os.path.isabs(manifest)))):
                    p.error("Invalid GWAS mode record: " + run_meta)
                if (requested_profile and requested_profile != profile) or (
                        requested_manifest and (not manifest or
                        os.path.abspath(requested_manifest) != manifest)):
                    p.error("Explicit GWAS mode conflicts with existing provenance; use a new output directory")
                if manifest and saved.get("manifest_sha256"):
                    if os.path.isfile(manifest):
                        with open(manifest, "rb") as f:
                            digest = hashlib.sha256(f.read()).hexdigest()
                        if digest != saved["manifest_sha256"]:
                            p.error("GWAS manifest content changed; use a new output directory")
                    elif requested_manifest:
                        p.error("Explicit GWAS manifest no longer exists")
                saved_build = saved.get("genome_build", "GRCh37")
                saved_policy = saved.get("filter_policy")
                saved_pop = saved.get("gwas_population") or ""
                if args.gwas_population is not None and saved_pop and args.gwas_population != saved_pop:
                    p.error("--gwas-population conflicts with recorded GWAS inputs; use a new output directory")
                if args.gwas_population is None and saved_pop:
                    args.gwas_population = saved_pop
                saved_skip = bool(saved.get("skip_meta", False))
                if args.skip_gwas_meta and not saved_skip:
                    p.error("--skip-gwas-meta conflicts with recorded meta-analysis mode; use a new output directory")
                if (saved_build not in ("GRCh37", "GRCh38")
                        or saved_policy not in ("none", "legacy")
                        or (profile and saved_policy != "legacy")):
                    p.error("Invalid GWAS identity/build/filter record: " + run_meta)
                for name, old, new in (("--genome-build", saved_build, args.genome_build),
                                       ("--gwas-filter-policy", saved_policy, args.gwas_filter_policy)):
                    if new is not None and new != old:
                        p.error(name + " conflicts with recorded GWAS inputs; use a new output directory")
                args.genome_build = saved_build
                args.gwas_filter_policy = saved_policy
                args.gwas_mode_recovered = True
                print("[gwas] recovered recorded mode from", run_meta)
            else:
                print("[gwas] recovered gwas_id={!r} from {}".format(saved_id, run_meta))
            args.gwas_id = saved_id
    if needs_mode and not profile and not manifest:
        p.error("GWAS steps require an explicit --gwas-profile (historical SLE "
                "reproduction, includes the known errors of the old "
                "implementation) or --gwas-manifest (corrected standard path); "
                "or reuse a successful recorded mode in " + run_meta)
    if args.gwas_id is None:
        if profile:
            args.gwas_id = "SLE"
        elif manifest:
            p.error("--gwas-manifest requires an explicit --gwas-id")
        elif needs_summary_id:
            p.error("summarizeTWAS requires --gwas-id when no recorded GWAS provenance exists in " + run_meta)
        else:
            args.gwas_id = ""
    args.genome_build = args.genome_build or "GRCh37"
    if profile:
        if args.genome_build != "GRCh37":
            p.error("Historical SLE profiles require --genome-build GRCh37")
        args.gwas_prefix = LEGACY_GWAS_PREFIX
        filter_policy = "legacy"
    elif manifest:
        is_single = args.skip_gwas_meta or (
            saved is not None and bool(saved.get("skip_meta", False))
        ) or ("alignGWAS" in args.steps and "metaGWAS" not in args.steps and "filterGWAS" not in args.steps)
        if is_single and (not args.gwas_mode_recovered or os.path.isfile(manifest)):
            import csv as _csv
            with open(manifest, "r", newline="") as _mf:
                _rows = list(_csv.DictReader(_mf, delimiter="\t"))
            if len(_rows) != 1:
                p.error("skipping metaGWAS/filterGWAS requires a single-study --gwas-manifest (found %d studies)" % len(_rows))
        args.gwas_prefix = args.gwas_id + ("_single" if is_single else "_meta")
        args.skip_meta_active = is_single
        filter_policy = "none"
    else:
        args.gwas_prefix = ""
        args.skip_meta_active = False
        filter_policy = ""
    if args.gwas_filter_policy is not None:
        if profile and args.gwas_filter_policy != "legacy":
            p.error("the legacy profile fixes --gwas-filter-policy legacy")
        filter_policy = args.gwas_filter_policy
    if needs_mode and saved and saved.get("prefix") != args.gwas_prefix:
        p.error("Recorded GWAS prefix disagrees with analysis identity")
    args.gwas_filter_policy = filter_policy
    args.metal_tbl = (os.path.join(gwas_dir, args.gwas_prefix + "1.tbl")
                      if args.gwas_prefix else "")
    args.gwas_filt = (os.path.abspath(args.gwas_summary) if args.gwas_summary
                      else (os.path.join(gwas_dir, args.gwas_prefix + "_filt.tsv")
                            if args.gwas_prefix else ""))
    args.ld_panel_id = args.ld_panel_id or args.gwas_population or "EUR_EAS"
    args.gwas_profile_active = profile
    args.gwas_manifest_active = manifest
    return args


def run(cmd, env):
    subprocess.run(cmd, env=env, check=True)


def conda_exe():
    """Resolve a conda new enough for modern conda-forge (old system conda breaks)."""
    exe = os.environ.get("CONDA_EXE")
    if exe:
        return exe
    miniforge = os.path.join(os.path.expanduser("~"), "miniforge3", "bin", "conda")
    if os.path.isfile(miniforge):
        return miniforge
    return "conda"


def with_bundled_r_libs(env_name, env):
    """Prepend bundled R packages for env_name (if staged under tools/R_libs) to R_LIBS_USER."""
    rlib = os.path.join(REPO, "tools", "R_libs", env_name)
    if os.path.isdir(rlib):
        env = env.copy()
        prev = env.get("R_LIBS_USER", "")
        env["R_LIBS_USER"] = rlib + (os.pathsep + prev if prev else "")
    return env


def resolve_lincs_env(env_base, requested):
    """Pick an R environment with fgsea installed unless explicitly provided."""
    if requested:
        return requested
    roots = []
    if env_base:
        roots.append(env_base)
    roots.append(os.path.join(os.path.dirname(conda_exe()), "..", "envs"))
    roots.append(os.path.join(os.path.expanduser("~"), "miniforge3", "envs"))
    for root in roots:
        r443_fgsea = os.path.join(root, "R4.4.3", "lib", "R", "library", "fgsea")
        if os.path.isdir(r443_fgsea):
            return "R4.4.3"
        r461_fgsea = os.path.join(root, "r_4.6.1", "lib", "R", "library", "fgsea")
        if os.path.isdir(r461_fgsea):
            return "r_4.6.1"
    return "R4.4.3"


def parse_lincs_args(argv):
    p = argparse.ArgumentParser(
        prog="medicell.py lincs",
        description="MediCell LINCS therapeutic prioritization: perturbational enrichment and cell-type scoring")
    p.add_argument("--output-dir", required=True,
                   help="output directory (may reuse a TWAS output-dir containing M01_deconvolution/M08_summary)")
    p.add_argument("--steps", nargs="+", metavar="STEP", default=None,
                   help="LINCS steps to run in order: " + ",".join(LINCS_STEPS)
                        + " (L01,L02,L03 aliases also work; default: all, or prepLincs,scoreLincs when --lincs-sig-csv is given)")
    for lstep in LINCS_STEPS:
        p.add_argument("--" + lstep, dest="flag_" + lstep, action="store_true",
                       help="run LINCS step " + lstep)
    p.add_argument("--trait", default="",
                   help="trait label used in output tags (default: auto from M04 gwas_run.json or 'trait')")
    p.add_argument("--tissue", "--label", dest="tissue", default="cohort",
                   help="cohort/compartment label used in output names (default cohort)")
    p.add_argument("--direction", choices=["pos", "neg", "both"], default="pos",
                   help="TWAS z-score direction to test: pos (>0), neg (<0), or both (default pos)")
    p.add_argument("--cell-types", default="",
                   help="optional space/comma-separated cell types to retain")

    twas = p.add_argument_group("TWAS and expression inputs")
    twas.add_argument("--twas-file", default="",
                      help="TWAS summary table (.csv/.tsv/.xlsx); defaults to M08_summary/allcelltype_result.csv")
    twas.add_argument("--twas-sheet", default="",
                      help="Excel sheet name when --twas-file is .xlsx")
    twas.add_argument("--twas-gene-col", default="",
                      help="explicit gene symbol column name in --twas-file")
    twas.add_argument("--twas-zscore-col", default="",
                      help="explicit z-score column name in --twas-file")
    twas.add_argument("--twas-padj-col", default="",
                      help="explicit adjusted p-value column name in --twas-file")
    twas.add_argument("--twas-celltype-col", default="",
                      help="explicit cell type column name in --twas-file")
    twas.add_argument("--twas-padj-cutoff", type=float, default=0.05,
                      help="TWAS adjusted p-value cutoff for fgsea gene selection (default 0.05)")
    twas.add_argument("--twas-zscore-cutoff", type=float, default=0.25,
                      help="minimum |z-score| cutoff for fgsea ranked genes (default 0.25)")
    twas.add_argument("--exp-matrix", default="",
                      help="pre-built gene x cell-type expression CSV (skips M01 aggregation)")
    twas.add_argument("--enigma-rdata", default="",
                      help="M01 ENIGMA RData path (default <output-dir>/M01_deconvolution/enigma.cse.trace_alpha_<alpha>.RData)")
    twas.add_argument("--alpha", type=float, default=0.1,
                      help="ENIGMA alpha used to locate default M01 RData (default 0.1)")
    twas.add_argument("--gene-pos", default="",
                      help="gene position CSV with gene_id,gene_name columns for mapping M01 rownames")

    lincs = p.add_argument_group("LINCS enrichment and scoring")
    lincs.add_argument("--lincs-gmt",
                       default=os.path.join(REPO, "data", "ref", "lincs", "LINCS_PRL_top250_human.gmt"),
                       help="LINCS/gene-set GMT file")
    lincs.add_argument("--lincs-pert-phase1",
                       default=os.path.join(REPO, "data", "ref", "lincs", "GSE92742_Broad_LINCS_pert_info.txt.gz"),
                       help="optional Phase 1 perturbation annotation (.txt.gz; pass '' to disable)")
    lincs.add_argument("--lincs-pert-phase2",
                       default=os.path.join(REPO, "data", "ref", "lincs", "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"),
                       help="optional Phase 2 perturbation annotation (.txt.gz; pass '' to disable)")
    lincs.add_argument("--lincs-input-csv", default="",
                       help="optional pre-split directional TWAS CSV for standalone L02")
    lincs.add_argument("--lincs-sig-csv", default="",
                       help="pre-computed LINCS enrichment CSV for L03 (skips L02 when --steps is omitted)")
    lincs.add_argument("--lincs-min-size", type=int, default=10,
                       help="fgsea minSize (default 10)")
    lincs.add_argument("--lincs-max-size", type=int, default=500,
                       help="fgsea maxSize (default 500)")
    lincs.add_argument("--lincs-sig-col", choices=["pval", "padj"], default="pval",
                       help="significance column used to filter LINCS hits (default pval)")
    lincs.add_argument("--lincs-sig-cutoff", type=float, default=0.05,
                       help="significance threshold for LINCS hits (default 0.05)")
    lincs.add_argument("--lincs-nes-filter", choices=["neg", "pos", "all"], default="neg",
                       help="NES direction filter for LINCS hits (default neg)")
    lincs.add_argument("--cores", type=int, default=None,
                       help="unified parallel cores/workers across LINCS steps (default 4)")
    lincs.add_argument("--lincs-cores", type=int, default=None,
                       help="alias/override for fgsea threads (defaults to --cores)")
    lincs.add_argument("--remove-ribo", action="store_true",
                       help="exclude ribosomal and mitochondrial genes before fgsea")
    lincs.add_argument("--top-n", type=int, default=50,
                       help="top LINCS perturbations retained per run (default 50)")
    lincs.add_argument("--lincs-rank-by", choices=["pval", "padj"], default="pval",
                       help="column used by L03 to rank perturbations (default pval)")
    lincs.add_argument("--exclude-pert-pattern", default=r"^BRD-|trt_sh\.css",
                       help="regex applied to 'pert_iname pert_type' to exclude perturbations (pass '' to disable)")

    env = p.add_argument_group("environment")
    env.add_argument("--env-base", default="",
                     help="portable environment root (unpacked conda-pack bundle)")
    env.add_argument("--lincs-env", default="",
                     help="conda env with fgsea/dplyr/tidyr/readxl (default: auto-detect)")

    args = p.parse_args(argv)

    if args.trait and not re.fullmatch(SAFE_NAME, args.trait):
        p.error("--trait must be a safe output label ([A-Za-z0-9_.-])")
    if not re.fullmatch(SAFE_NAME, args.tissue):
        p.error("--tissue/--label must be a safe output label ([A-Za-z0-9_.-])")

    lincs_flagged = [s for s in LINCS_STEPS if getattr(args, "flag_" + s)]
    if lincs_flagged and args.steps is not None:
        p.error("step flags cannot be combined with --steps")
    if lincs_flagged:
        args.steps = lincs_flagged
    elif args.steps is None:
        if args.lincs_sig_csv:
            args.steps = ["prepLincs", "scoreLincs"]
        else:
            args.steps = LINCS_STEPS.copy()
    else:
        selected = []
        for value in args.steps:
            for part in value.split(","):
                if not part.strip():
                    p.error("--steps cannot contain empty values; choose from: " + ",".join(LINCS_STEPS))
                selected.extend(part.split())
        canonical = []
        unknown = []
        for s in selected:
            c = LINCS_STEP_ALIASES.get(s, s)
            if c not in LINCS_STEPS:
                unknown.append(s)
            elif c not in canonical:
                canonical.append(c)
        if unknown:
            p.error("unknown LINCS steps: " + ",".join(sorted(set(unknown)))
                    + "; choose from: " + ",".join(LINCS_STEPS))
        args.steps = [s for s in LINCS_STEPS if s in canonical]

    cell_types = []
    if args.cell_types:
        cell_types = args.cell_types.replace(",", " ").split()
        if any(not re.fullmatch(SAFE_NAME, ct) for ct in cell_types):
            p.error("cell type names must be safe output labels ([A-Za-z0-9_.-])")
        if len(set(cell_types)) != len(cell_types):
            p.error("--cell-types contains duplicates")
    args.cell_types = cell_types

    if "scoreLincs" in args.steps and len(args.cell_types) == 1:
        p.error("LINCS scoring (scoreLincs) requires at least 2 cell types; single-cell-type input yields undefined cross-cell-type expression z-scores")
    if args.top_n < 1:
        p.error("--top-n must be >= 1")
    if args.lincs_min_size < 1 or args.lincs_max_size < args.lincs_min_size:
        p.error("--lincs-min-size and --lincs-max-size must satisfy 1 <= min <= max")
    if args.cores is not None and args.cores < 1:
        p.error("--cores must be >= 1")
    if args.lincs_cores is not None and args.lincs_cores < 1:
        p.error("--lincs-cores must be >= 1")
    args.cores_resolved = (args.cores if args.cores is not None
                           else (args.lincs_cores if args.lincs_cores is not None else 4))
    args.lincs_cores_resolved = (args.lincs_cores if args.lincs_cores is not None
                                 else args.cores_resolved)
    if not (0 < args.twas_padj_cutoff <= 1) or not (0 < args.lincs_sig_cutoff <= 1):
        p.error("p-value cutoffs must be in (0, 1]")
    if args.twas_zscore_cutoff < 0 or not math.isfinite(args.twas_zscore_cutoff):
        p.error("--twas-zscore-cutoff must be finite and >= 0")

    return args


def run_lincs_mode(argv):
    args = parse_lincs_args(argv)
    env_base = os.path.abspath(args.env_base) if args.env_base else ""
    lincs_env = resolve_lincs_env(env_base, args.lincs_env)
    print_lock = threading.Lock()

    def conda_run(env_name, cmd, env):
        env = with_bundled_r_libs(env_name, env)
        if env_base:
            env = env.copy()
            bindir = os.path.join(env_base, env_name, "bin")
            env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
            with print_lock:
                print("[run]", env_name, " ".join(cmd), flush=True)
            subprocess.run(cmd, env=env, check=True)
        else:
            full = [conda_exe(), "run", "-n", env_name, "--no-capture-output"] + cmd
            with print_lock:
                print("[run]", " ".join(full), flush=True)
            subprocess.run(full, env=env, check=True)

    out = os.path.abspath(args.output_dir)
    dirs = {
        "l01": os.path.join(out, "L01_inputs"),
        "l02": os.path.join(out, "L02_enrichment"),
        "l03": os.path.join(out, "L03_drug_scores"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    def abspath_opt(val):
        return os.path.abspath(val) if val else ""

    trait = args.trait
    if not trait:
        gwas_meta = read_metadata(os.path.join(out, "M04_gwas", "gwas_run.json"))
        if isinstance(gwas_meta, dict) and isinstance(gwas_meta.get("gwas_id"), str) and gwas_meta["gwas_id"]:
            trait = gwas_meta["gwas_id"]
        else:
            trait = "trait"

    twas_file = abspath_opt(args.twas_file)
    if not twas_file:
        m08_all = os.path.join(out, "M08_summary", "allcelltype_result.csv")
        m08_fdr = os.path.join(out, "M08_summary", "FDRsig_allcelltype_result.csv")
        if os.path.isfile(m08_all):
            twas_file = m08_all
        elif os.path.isfile(m08_fdr):
            twas_file = m08_fdr

    enigma_rdata = abspath_opt(args.enigma_rdata)
    if not enigma_rdata and not args.exp_matrix:
        default_m01 = os.path.join(out, "M01_deconvolution",
                                   "enigma.cse.trace_alpha_{}.RData".format(args.alpha))
        if os.path.isfile(default_m01):
            enigma_rdata = default_m01

    if "prepLincs" in args.steps:
        if not twas_file:
            sys.exit("L01 requires --twas-file (or M08_summary/allcelltype_result.csv in --output-dir)")
        if not args.exp_matrix and not enigma_rdata:
            sys.exit("L01 requires --exp-matrix or --enigma-rdata (or M01_deconvolution RData in --output-dir)")

    if "enrichLincs" in args.steps:
        if not os.path.isfile(args.lincs_gmt):
            sys.exit("L02 requires existing file for --lincs-gmt: {}".format(args.lincs_gmt))
        for flag, fpath in (("--lincs-pert-phase1", args.lincs_pert_phase1),
                            ("--lincs-pert-phase2", args.lincs_pert_phase2)):
            if fpath and not os.path.isfile(fpath):
                sys.exit("L02 annotation file not found for {}: {}".format(flag, fpath))

    directions = ["pos", "neg"] if args.direction == "both" else [args.direction]
    if args.direction == "both" and (args.lincs_sig_csv or args.lincs_input_csv):
        sys.exit("--lincs-sig-csv / --lincs-input-csv cannot be combined with --direction both")

    def lincs_script(step):
        return os.path.join(SCRIPTS, LINCS_SCRIPT_FILES[step])

    total_cores = args.cores_resolved
    lincs_cores = args.lincs_cores_resolved
    dir_workers = min(len(directions), total_cores)
    l02_inner_cores = max(1, lincs_cores // dir_workers)

    print("[lincs steps]", " -> ".join(args.steps))
    print("[lincs trait/label]", trait, "/", args.tissue, "directions:", ",".join(directions), "cores:", total_cores)
    print("[lincs env]", lincs_env)

    def make_env(direction, step_cores):
        twas_tag = "{}_{}_{}".format(trait, args.tissue, direction)
        l01_dir_csv = os.path.join(dirs["l01"], "{}.csv".format(twas_tag))
        l01_twas_tsv = os.path.join(dirs["l01"], "{}_{}_twas_all.tsv".format(trait, args.tissue))
        l01_exp_csv = os.path.join(dirs["l01"], "{}_exp_matrix.csv".format(args.tissue))
        l02_sig_csv = os.path.join(dirs["l02"], "LINCS_sig_{}.csv".format(twas_tag))

        l02_in_csv = abspath_opt(args.lincs_input_csv) if args.lincs_input_csv else l01_dir_csv
        l03_sig_csv = abspath_opt(args.lincs_sig_csv) if args.lincs_sig_csv else l02_sig_csv
        l03_twas_tsv = l01_twas_tsv if ("prepLincs" in args.steps or os.path.isfile(l01_twas_tsv)) else twas_file
        l03_exp_csv = l01_exp_csv if ("prepLincs" in args.steps or os.path.isfile(l01_exp_csv)) else abspath_opt(args.exp_matrix)

        if "scoreLincs" in args.steps and (not l03_twas_tsv or not l03_exp_csv):
            sys.exit("L03 requires L01 outputs or explicit --twas-file and --exp-matrix")

        env = os.environ.copy()
        env.update({
            "MC_CORES": str(step_cores),
            "MC_LINCS_L01_DIR": dirs["l01"],
            "MC_LINCS_L02_DIR": dirs["l02"],
            "MC_LINCS_L03_DIR": dirs["l03"],
            "MC_LINCS_TRAIT": trait,
            "MC_LINCS_LABEL": args.tissue,
            "MC_LINCS_DIRECTION": direction,
            "MC_LINCS_TWAS_LABEL": twas_tag,
            "MC_CELLTYPES": " ".join(args.cell_types),
            "MC_LINCS_MIN_CTS": "2" if "scoreLincs" in args.steps else "1",
            "MC_LINCS_TWAS_FILE": twas_file,
            "MC_LINCS_TWAS_SHEET": args.twas_sheet,
            "MC_LINCS_GENE_COL": args.twas_gene_col,
            "MC_LINCS_Z_COL": args.twas_zscore_col,
            "MC_LINCS_PADJ_COL": args.twas_padj_col,
            "MC_LINCS_CT_COL": args.twas_celltype_col,
            "MC_LINCS_L03_GENE_COL": "" if l03_twas_tsv == l01_twas_tsv else args.twas_gene_col,
            "MC_LINCS_L03_CT_COL": "" if l03_twas_tsv == l01_twas_tsv else args.twas_celltype_col,
            "MC_LINCS_L03_Z_COL": "" if l03_twas_tsv == l01_twas_tsv else args.twas_zscore_col,
            "MC_LINCS_TWAS_PADJ_CUTOFF": str(args.twas_padj_cutoff),
            "MC_LINCS_ZSCORE_CUTOFF": str(args.twas_zscore_cutoff),
            "MC_LINCS_EXP_MATRIX": abspath_opt(args.exp_matrix),
            "MC_ENIGMA_RDATA": enigma_rdata,
            "MC_GENE_POS": abspath_opt(args.gene_pos),
            "MC_LINCS_DIR_CSV": l02_in_csv,
            "MC_LINCS_GMT": abspath_opt(args.lincs_gmt),
            "MC_LINCS_PERT_PHASE1": abspath_opt(args.lincs_pert_phase1),
            "MC_LINCS_PERT_PHASE2": abspath_opt(args.lincs_pert_phase2),
            "MC_LINCS_MIN_SIZE": str(args.lincs_min_size),
            "MC_LINCS_MAX_SIZE": str(args.lincs_max_size),
            "MC_LINCS_SIG_COL": args.lincs_sig_col,
            "MC_LINCS_SIG_CUTOFF": str(args.lincs_sig_cutoff),
            "MC_LINCS_NES_FILTER": args.lincs_nes_filter,
            "MC_LINCS_CORES": str(step_cores),
            "MC_LINCS_REMOVE_RIBO": "1" if args.remove_ribo else "0",
            "MC_LINCS_SIG_CSV": l03_sig_csv,
            "MC_LINCS_TWAS_TSV": l03_twas_tsv,
            "MC_LINCS_EXP_CSV": l03_exp_csv,
            "MC_LINCS_TOP_N": str(args.top_n),
            "MC_LINCS_RANK_BY": args.lincs_rank_by,
            "MC_LINCS_EXCLUDE_PATTERN": args.exclude_pert_pattern if args.exclude_pert_pattern else "NONE",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        })
        return env

    if "prepLincs" in args.steps:
        for direction in directions:
            conda_run(lincs_env, ["Rscript", lincs_script("prepLincs")], make_env(direction, total_cores))

    def run_dir_step(step, cmd_runner, per_dir_cores):
        if dir_workers <= 1:
            for direction in directions:
                conda_run(lincs_env, cmd_runner(step), make_env(direction, per_dir_cores))
        else:
            with ThreadPoolExecutor(max_workers=dir_workers) as pool:
                futs = [pool.submit(conda_run, lincs_env, cmd_runner(step), make_env(d, per_dir_cores))
                        for d in directions]
                for fut in futs:
                    fut.result()

    if "enrichLincs" in args.steps:
        run_dir_step("enrichLincs", lambda s: ["Rscript", lincs_script(s)], l02_inner_cores)
    if "scoreLincs" in args.steps:
        run_dir_step("scoreLincs", lambda s: ["bash", lincs_script(s)], 1)

    print("MediCell LINCS run complete. Results in:", out)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "lincs":
        run_lincs_mode(sys.argv[2:])
        return
    if "--lincs" in sys.argv[1:]:
        run_lincs_mode([a for a in sys.argv[1:] if a != "--lincs"])
        return
    argv = sys.argv[1:]
    if argv and argv[0] == "cell-type-twas":
        argv = ["--cell-type-twas"] + argv[1:]
    elif argv and (argv[0] in ALL_STEPS or argv[0] in STEP_ALIASES):
        argv = ["--steps", argv[0]] + argv[1:]
    args = parse_args(argv)
    env_base = os.path.abspath(args.env_base) if args.env_base else ""
    metal_default = os.path.join(REPO, "tools", "METAL", "metal")
    if env_base and args.metal_bin == metal_default:
        args.metal_bin = os.path.join(env_base, "r3.6.3", "bin", "metal")

    print_lock = threading.Lock()

    def conda_run(env_name, cmd, env):
        env = with_bundled_r_libs(env_name, env)
        if env_base:
            # portable mode: run binaries directly, env bin first on PATH
            env = env.copy()
            bindir = os.path.join(env_base, env_name, "bin")
            env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
            with print_lock:
                print("[run]", env_name, " ".join(cmd), flush=True)
            subprocess.run(cmd, env=env, check=True)
        else:
            full = [conda_exe(), "run", "-n", env_name, "--no-capture-output"] + cmd
            with print_lock:
                print("[run]", " ".join(full), flush=True)
            subprocess.run(full, env=env, check=True)

    out = os.path.abspath(args.output_dir)
    chrs = args.chromosomes
    cell_types = args.cell_types
    peer_nk = args.peer_nk_resolved
    peer_keep = args.peer_keep_resolved
    steps = args.steps
    if args.step_aliases:
        print("[step aliases]",
              ", ".join("%s->%s" % pair for pair in args.step_aliases))
    print("[steps]", " -> ".join(steps))
    if args.gwas_profile_active:
        print("[gwas] WARNING: legacy profile %r reproduces the historical SLE "
              "implementation, including its known statistical errors "
              "(mis-projected beta/SE columns, swapped effect-allele declaration); "
              "use --gwas-manifest for the corrected standard path"
              % args.gwas_profile_active)

    dirs = {
        "m01": os.path.join(out, "M01_deconvolution"),
        "exp": os.path.join(out, "M02_expression"),
        "cov": os.path.join(out, "M02_covariates"),
        "m03": os.path.join(out, "M03_models"),
        "db": os.path.join(out, "M06_databases"),
        "gwas": os.path.join(out, "M04_gwas"),
        "ld_dosage": (args.ld_dosage_dir if args.ld_dosage_dir
                      else os.path.join(out, "M05_ld", "dosage")),
        "ld_cov": os.path.join(out, "M05_ld", "cov_" + args.ld_panel_id),
        "twas": os.path.join(out, "M07_twas"),
        "summary": os.path.join(out, "M08_summary"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    def record_gwas_mode():
        if args.gwas_mode_recovered:
            return
        digest = ""
        if args.gwas_manifest_active and os.path.isfile(args.gwas_manifest_active):
            with open(args.gwas_manifest_active,"rb") as f:
                digest = hashlib.sha256(f.read()).hexdigest()
        write_metadata(os.path.join(dirs["gwas"], "gwas_run.json"), {
            "manifest_sha256": digest,
            "profile": args.gwas_profile_active or "",
            "manifest": (os.path.abspath(args.gwas_manifest_active)
                         if args.gwas_manifest_active else ""),
            "gwas_id": args.gwas_id,
            "genome_build": args.genome_build,
            "prefix": args.gwas_prefix,
            "metal_tbl": args.metal_tbl,
            "filtered": args.gwas_filt,
            "filter_policy": args.gwas_filter_policy,
            "gwas_population": args.gwas_population or "",
            "skip_meta": bool(getattr(args, "skip_meta_active", False)),
        })

    enigma_rdata = os.path.join(dirs["m01"],
                                "enigma.cse.trace_alpha_{}.RData".format(args.alpha))

    def input_path(value):
        return os.path.abspath(value) if value else ""

    base = {
        "MC_ALPHA": str(args.alpha),
        "MC_TPM_MIN": str(args.tpm_min),
        "MC_MIN_PCT": str(args.min_pct),
        "MC_MIN_SAMPLES": str(args.min_samples),
        "MC_BULK_TPM": input_path(args.bulk_tpm),
        "MC_ENIGMA_REF": input_path(args.enigma_ref),
        "MC_CBS_FRAX": input_path(args.cibersortx_frax),
        "MC_COVARIATES": input_path(args.covariates),
        "MC_CAT_COVS": args.cat_covs,
        "MC_MIN_LEVEL_COUNT": str(args.min_level_count),
        "MC_PEER_NK": "" if peer_nk is None else str(peer_nk),
        "MC_PEER_KEEP": "" if peer_keep is None else str(peer_keep),
        "MC_SNP_ANNOT_PATTERN": input_path(args.snp_annot_pattern),
        "MC_GENOTYPE_PATTERN": input_path(args.genotype_pattern),
        "MC_GENE_ANNOT": input_path(args.gene_annot),
        "MC_GENE_ANNOT_FORMAT": args.gene_annot_format,
        "MC_SNP_ANNOT": input_path(args.snp_annot_file),
        "MC_GWAS_IN_DIR": input_path(args.gwas_dir),
        "MC_GWAS_PROFILE": args.gwas_profile_active or "",
        "MC_GWAS_MANIFEST": (os.path.abspath(args.gwas_manifest_active)
                             if args.gwas_manifest_active else ""),
        "MC_GWAS_ID": args.gwas_id,
        "MC_GWAS_PREFIX": args.gwas_prefix,
        "MC_GWAS_OUT": dirs["gwas"],
        "MC_GWAS_FILTER_POLICY": args.gwas_filter_policy,
        "MC_METAL_TBL": args.metal_tbl,
        "MC_GWAS_FILT": args.gwas_filt,
        "MC_GWAS_SINGLE_OUT": (args.gwas_filt if getattr(args, "skip_meta_active", False) else ""),
        "MC_GWAS_POPULATION": args.gwas_population or "",
        "MC_G1000_VCF_PATTERN": input_path(args.g1000_vcf_pattern),
        "MC_G1000_SAMPLES": input_path(args.g1000_samples),
        "MC_GENE_POS": input_path(args.gene_pos),
        "MC_SAMPLE_ID_POLICY": args.sample_id_policy,
        "MC_COV_LAYOUT": args.cov_layout,
        "MC_GENOME_BUILD": args.genome_build,
        "MC_TRAIN_POPULATION": args.training_population,
        "MC_MAF": str(args.maf),
        "MC_MAF_POLICY": args.maf_policy,
        "MC_LD_PANEL_ID": args.ld_panel_id,
        "MC_TWAS_MODEL_SET": args.twas_model_set or "",
        "MC_CORES": str(args.cores_resolved),
        "MC_SEED": "" if args.seed is None else str(args.seed),
        "MC_TRAINING_CORES": str(1 if args.seed is not None
                                 else (args.training_cores or args.cores_resolved)),
        "MC_METAXCAN_DIR": os.path.abspath(args.metaxcan_dir),
        "MC_METAL_BIN": os.path.abspath(args.metal_bin),
        "MC_BCFTOOLS": args.bcftools,
        "MC_BCF_PLUGINS": args.bcf_plugins,
        "MC_ELNET_SCRIPT": os.path.join(SCRIPTS, "helpers", "gtex_v7_nested_cv_elnet.R"),
        "MC_CELLTYPES": " ".join(cell_types),
        "MC_CHROMS": " ".join(str(c) for c in chrs),
        "MC_CHROMS_CSV": ",".join(str(c) for c in chrs),
        "MC_LD_SKIP_EXTRACT": "1" if args.ld_dosage_dir else "0",
    }

    def mc(**overrides):
        env = os.environ.copy()
        env.update(base)
        env.update(overrides)
        return env

    def script(step):
        return os.path.join(SCRIPTS, SCRIPT_FILES[step])

    total_cores = args.cores_resolved
    blas_single = {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    blas_multi = {
        "OMP_NUM_THREADS": str(total_cores),
        "OPENBLAS_NUM_THREADS": str(total_cores),
        "MKL_NUM_THREADS": str(total_cores),
    }

    if "inferExp" in steps:
        conda_run(args.r_env, ["Rscript", script("inferExp")],
                  mc(MC_OUT_DIR=dirs["m01"], **blas_multi))

    if "normExp" in steps:
        conda_run(args.peer_env, ["Rscript", script("normExp")],
                  mc(MC_ENIGMA_RDATA=enigma_rdata, MC_EXP_DIR=dirs["exp"], **blas_single))

    if "peerCov" in steps:
        conda_run(args.peer_env, ["Rscript", script("peerCov")],
                  mc(MC_EXP_DIR=dirs["exp"], MC_COV_DIR=dirs["cov"], **blas_single))

    if "trainModel" in steps:
        m03_jobs = [(ct, chr_) for ct in cell_types for chr_ in chrs]
        if args.seed is not None:
            m03_inner = 1
            m03_workers = min(len(m03_jobs), total_cores)
        elif args.training_cores is not None:
            m03_inner = args.training_cores
            m03_workers = max(1, min(len(m03_jobs), total_cores // m03_inner))
        else:
            m03_workers = min(len(m03_jobs), total_cores)
            m03_inner = max(1, total_cores // m03_workers)

        def run_m03_one(ct, chr_):
            conda_run(args.r_env, ["Rscript", script("trainModel"), ct, str(chr_)],
                      mc(MC_COV_FILE=os.path.join(dirs["cov"], ct + "_covariates.txt"),
                         MC_EXP_FILE=os.path.join(dirs["exp"], ct + "_transformed_expression.txt"),
                         MC_M03_OUT=dirs["m03"],
                         MC_TRAINING_CORES=str(m03_inner),
                         **blas_single))

        if m03_workers <= 1:
            for ct, chr_ in m03_jobs:
                run_m03_one(ct, chr_)
        else:
            with ThreadPoolExecutor(max_workers=m03_workers) as pool:
                futs = [pool.submit(run_m03_one, ct, chr_) for ct, chr_ in m03_jobs]
                for fut in futs:
                    fut.result()

    if "buildDB" in steps:
        conda_run(args.r_env, ["Rscript", script("buildDB")],
                  mc(MC_M03_SUMMARY_DIR=os.path.join(dirs["m03"], "summary"),
                     MC_M03_WEIGHTS_DIR=os.path.join(dirs["m03"], "weights"),
                     MC_DB_DIR=dirs["db"],
                     MC_CHROMS=base["MC_CHROMS_CSV"],
                     **blas_single))

    if "alignGWAS" in steps:
        conda_run(args.r_env, ["Rscript", script("alignGWAS")], mc())
        record_gwas_mode()

    if "metaGWAS" in steps:
        conda_run(args.r_env, ["bash", script("metaGWAS")],
                  mc(MC_GWAS_OUT=dirs["gwas"], MC_METAL_TBL=args.metal_tbl))
        record_gwas_mode()

    if "filterGWAS" in steps:
        conda_run(args.r_env, ["Rscript", script("filterGWAS")],
                  mc(MC_METAL_TBL=args.metal_tbl,
                     MC_GWAS_OUT=dirs["gwas"],
                     MC_GWAS_FILT=args.gwas_filt))
        record_gwas_mode()

    if "extractLD" in steps:
        conda_run(args.py_env, ["bash", script("extractLD")],
                  mc(MC_LD_DOSAGE_DIR=dirs["ld_dosage"],
                     MC_LD_COV_DIR=dirs["ld_cov"],
                     **blas_single))

    if "buildLD" in steps:
        conda_run(args.py_env, ["bash", script("buildLD")],
                  mc(MC_DB_DIR=dirs["db"], MC_LD_COV_DIR=dirs["ld_cov"],
                     **blas_single))

    if "runTWAS" in steps:
        conda_run(args.py_env, ["bash", script("runTWAS")],
                  mc(MC_DB_DIR=dirs["db"], MC_LD_COV_DIR=dirs["ld_cov"],
                     MC_GWAS_FILT=args.gwas_filt, MC_TWAS_OUT=dirs["twas"],
                     **blas_single))

    if "summarizeTWAS" in steps:
        conda_run(args.r_env, ["Rscript", script("summarizeTWAS")],
                  mc(MC_TWAS_DIR=dirs["twas"], MC_SUMMARY_OUT=dirs["summary"],
                     MC_GWAS_ID=args.gwas_id, **blas_single))

    print("MediCell run complete. Results in:", out)


if __name__ == "__main__":
    main()
