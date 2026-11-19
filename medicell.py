#!/usr/bin/env python3
"""MediCell pipeline driver for cell-type TWAS and LINCS therapeutic prioritization.

Modes:
  1. Cell-type TWAS (`medicell.py --cell-type-twas` / `medicell.py cell-type-twas`)
  2. LINCS prioritization (`medicell.py lincs` / `medicell.py --lincs`)

All step paths and parameters are passed to child scripts via MC_* environment
variables, with outputs organized under --output-dir.

Canonical TWAS step order:
    inferExp -> normExp -> peerCov -> trainModel -> buildDB ->
    alignGWAS -> metaGWAS -> filterGWAS -> extractLD -> buildLD ->
    runTWAS -> summarizeTWAS
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import threading

REPO = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(REPO, "script")

# ---------------------------------------------------------------------------
# Step registries and constants
# ---------------------------------------------------------------------------

TWAS_STEPS = (
    ("inferExp", "M01_ENIGMA.R", "deconvolve bulk TPM into cell-type expression (ENIGMA)"),
    ("normExp", "M02_DataProcess_exp_fordb.R", "filter expressed genes and INT-transform per cell type"),
    ("peerCov", "M02_DataProcess_covariance_for_db.R", "estimate PEER hidden factors and format covariates"),
    ("trainModel", "M03_ModelTraining.R", "train cell-type elastic-net GReX models"),
    ("buildDB", "M04_Makedb.R", "build PrediXcan SQLite weight databases"),
    ("alignGWAS", "M05_DataProcess_GWAS_align4METAL.R", "align GWAS summary statistics to the SNP panel"),
    ("metaGWAS", "M05_DataProcess_GWAS_metal.sh", "run METAL GWAS meta-analysis"),
    ("filterGWAS", "M05_DataProcess_GWAS_metal_filter.R", "filter meta-analyzed GWAS variants"),
    ("extractLD", "M06_LD.sh", "extract and align 1000G reference dosages"),
    ("buildLD", "M06_Makeld.sh", "compute cell-type LD covariance matrices"),
    ("runTWAS", "M07_TWAS.sh", "run S-PrediXcan cell-type TWAS"),
    ("summarizeTWAS", "M08_ResultStatistic.R", "summarize TWAS results with FDR and Bonferroni"),
)
ALL_STEPS = [name for name, _, _ in TWAS_STEPS]
STEP_DESCRIPTIONS = {name: desc for name, _, desc in TWAS_STEPS}
SCRIPT_FILES = {name: script for name, script, _ in TWAS_STEPS}

# Retained for backward-compatible --steps CLI invocation.
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

LINCS_STEP_SPECS = (
    ("prepLincs", "L01_lincs_prepare.R"),
    ("enrichLincs", "L02_lincs_enrichment.R"),
    ("scoreLincs", "L03_lincs_drug_score.sh"),
)
LINCS_STEPS = [name for name, _ in LINCS_STEP_SPECS]
LINCS_SCRIPT_FILES = dict(LINCS_STEP_SPECS)
LINCS_STEP_ALIASES = {"L01": "prepLincs", "L02": "enrichLincs", "L03": "scoreLincs"}

GWAS_STEPS = {"alignGWAS", "metaGWAS", "filterGWAS"}
GWAS_SUMMARY_REQUIRED_COLS = {"MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value"}
LEGACY_GWAS_PREFIX = "SLE_metal_096_917_183"
LEGACY_GWAS_PROFILES = ("sle-legacy-pbmc", "sle-legacy-kidney")
SAFE_NAME = r"[A-Za-z0-9][A-Za-z0-9_.\-]*"


# ---------------------------------------------------------------------------
# Shared runtime and I/O helpers
# ---------------------------------------------------------------------------

def abspath_opt(path):
    """Return absolute path when non-empty, else empty string."""
    return os.path.abspath(path) if path else ""


def ensure_dirs(mapping):
    """Create every directory in mapping and return mapping."""
    for d in mapping.values():
        os.makedirs(d, exist_ok=True)
    return mapping


def sha256_file(path):
    """Compute SHA-256 hex digest of a file."""
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def write_metadata(path, data):
    """Atomically write JSON metadata next to its destination path."""
    with tempfile.NamedTemporaryFile(
        mode="w", dir=os.path.dirname(path), prefix=".gwas-run-", delete=False
    ) as f:
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
    """Read JSON metadata if present and valid, else return None."""
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def conda_exe():
    """Locate conda executable ($CONDA_EXE, ~/miniforge3/bin/conda, or PATH)."""
    exe = os.environ.get("CONDA_EXE")
    if exe:
        return exe
    miniforge = os.path.join(os.path.expanduser("~"), "miniforge3", "bin", "conda")
    return miniforge if os.path.isfile(miniforge) else "conda"


def with_bundled_r_libs(env_name, env):
    """Prepend staged tools/R_libs/<env_name> to R_LIBS_USER when available."""
    rlib = os.path.join(REPO, "tools", "R_libs", env_name)
    if os.path.isdir(rlib):
        env = env.copy()
        prev = env.get("R_LIBS_USER", "")
        env["R_LIBS_USER"] = rlib + (os.pathsep + prev if prev else "")
    return env


def make_conda_runner(env_base):
    """Build a thread-safe subprocess runner supporting conda or conda-pack env_base."""
    print_lock = threading.Lock()

    def conda_run(env_name, cmd, env):
        env = with_bundled_r_libs(env_name, env)
        if env_base:
            env = env.copy()
            bindir = os.path.join(env_base, env_name, "bin")
            env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
            full = cmd
            label = [env_name] + cmd
        else:
            full = [conda_exe(), "run", "-n", env_name, "--no-capture-output"] + cmd
            label = full
        with print_lock:
            print("[run]", " ".join(label), flush=True)
        subprocess.run(full, env=env, check=True)

    return conda_run


def run_parallel(workers, items, fn):
    """Execute fn(item) serially or in a ThreadPoolExecutor depending on worker count."""
    if workers <= 1:
        for item in items:
            fn(item)
        return
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for fut in [pool.submit(fn, item) for item in items]:
            fut.result()


def script_cmd(path, *extra_args):
    """Build command list using bash for .sh scripts and Rscript for .R scripts."""
    runner = "bash" if path.endswith(".sh") else "Rscript"
    return [runner, path] + [str(a) for a in extra_args]


def blas_env(threads):
    """Return single- or multi-threaded BLAS environment settings."""
    t = str(threads)
    return {"OMP_NUM_THREADS": t, "OPENBLAS_NUM_THREADS": t, "MKL_NUM_THREADS": t}


def resolve_lincs_env(env_base, requested):
    """Pick an R environment with fgsea installed unless explicitly provided."""
    if requested:
        return requested
    roots = [env_base] if env_base else []
    roots.append(os.path.join(os.path.dirname(conda_exe()), "..", "envs"))
    roots.append(os.path.join(os.path.expanduser("~"), "miniforge3", "envs"))
    for root in roots:
        for env_name in ("R4.4.3", "r_4.6.1"):
            if os.path.isdir(os.path.join(root, env_name, "lib", "R", "library", "fgsea")):
                return env_name
    return "R4.4.3"


# ---------------------------------------------------------------------------
# Shared CLI validation helpers
# ---------------------------------------------------------------------------

def validate_safe_label(parser, flag, value):
    """Ensure value matches SAFE_NAME when non-empty."""
    if value and not re.fullmatch(SAFE_NAME, value):
        parser.error(f"{flag} must be a safe output label ([A-Za-z0-9_.-])")


def check_finite_bounds(parser, value, flag, *, gt=None, ge=None, le=None):
    """Validate that a numeric option is finite and within specified bounds."""
    if not math.isfinite(value):
        parser.error(f"{flag} must be finite")
    if gt is not None and value <= gt:
        parser.error(f"{flag} must be > {gt}")
    if ge is not None and value < ge:
        parser.error(f"{flag} must be >= {ge}")
    if le is not None and value > le:
        parser.error(f"{flag} must be <= {le}")


def parse_cell_types(parser, raw, required=True):
    """Parse comma/space-separated cell types and enforce safe, unique names."""
    cell_types = raw.replace(",", " ").split() if raw else []
    if not cell_types:
        if required:
            parser.error("--cell-types cannot be empty")
        return []
    if any(not re.fullmatch(SAFE_NAME, ct) for ct in cell_types):
        parser.error("cell type names must be safe output labels ([A-Za-z0-9_.-])")
    if len(set(cell_types)) != len(cell_types):
        parser.error("--cell-types contains duplicates")
    return cell_types


def parse_chromosomes(parser, raw):
    """Parse comma-separated autosome list (1-22), preserving order without duplicates."""
    chrs = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            parser.error("--chromosomes cannot contain empty values")
        if not re.fullmatch(r"\d+", token) or not 1 <= int(token) <= 22:
            parser.error(f"invalid chromosome {token!r} (must be 1-22)")
        num = int(token)
        if num not in chrs:
            chrs.append(num)
    return chrs


def resolve_step_selection(
    parser, args, canonical_steps, aliases, default_steps_fn,
    *, empty_msg, unknown_prefix, has_full_flag=False
):
    """Normalize --steps / --<step> flags into ordered canonical step names."""
    flagged = [step for step in canonical_steps if getattr(args, "flag_" + step)]
    full_active = has_full_flag and args.cell_type_twas
    if flagged and (full_active or args.steps is not None):
        suffix = "--cell-type-twas or --steps" if has_full_flag else "--steps"
        parser.error(f"step shortcut flags (--{',--'.join(flagged)}) cannot be combined with {suffix}")

    if flagged:
        return flagged, []
    if args.steps is None:
        return default_steps_fn(args), []

    selected = []
    for value in args.steps:
        for part in value.split(","):
            if not part.strip():
                parser.error(f"{empty_msg}; choose from: {','.join(canonical_steps)}")
            selected.extend(part.split())

    canonical, unknown, used_aliases = [], [], []
    for s in selected:
        c = aliases.get(s, s)
        if c not in canonical_steps:
            unknown.append(s)
            continue
        if s != c and (s, c) not in used_aliases:
            used_aliases.append((s, c))
        if c not in canonical:
            canonical.append(c)

    if unknown:
        parser.error(f"{unknown_prefix}: {','.join(sorted(set(unknown)))}; choose from: {','.join(canonical_steps)}")
    return [step for step in canonical_steps if step in canonical], used_aliases


def resolve_cores(parser, primary, secondary, secondary_flag, default=4):
    """Validate and resolve --cores and an optional secondary step-core override."""
    if primary is not None and primary < 1:
        parser.error("--cores must be >= 1")
    if secondary is not None and secondary < 1:
        parser.error(f"{secondary_flag} must be >= 1")
    return primary if primary is not None else (secondary if secondary is not None else default)


# ---------------------------------------------------------------------------
# Cell-type TWAS argument parsing and GWAS provenance
# ---------------------------------------------------------------------------

def build_twas_parser():
    """Construct the argparse parser for the cell-type TWAS workflow."""
    p = argparse.ArgumentParser(
        description="MediCell pipeline driver (use `medicell.py lincs --help` for LINCS therapeutic prioritization)"
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--alpha", type=float, default=0.1, help="ENIGMA alpha (default 0.1)")
    p.add_argument("--tpm-min", type=float, default=0.1, help="normExp gene filter: min inferred TPM (default 0.1)")
    p.add_argument("--min-pct", type=float, default=0.01, help="normExp gene filter: min sample fraction, 0-1 (default 0.01)")
    p.add_argument("--min-samples", type=int, default=5, help="normExp gene filter: min sample count floor (default 5)")
    p.add_argument("--cat-covs", default="", help="peerCov: comma-separated categorical covariate columns to dummy-encode")
    p.add_argument("--min-level-count", type=int, default=5, help="peerCov: drop categorical levels with fewer samples (default 5)")
    p.add_argument("--peer-nk", type=int, default=None, help="PEER hidden factors; required explicitly for peerCov")
    p.add_argument("--peer-keep", type=int, default=None, help="PEER factors kept; required explicitly for peerCov")

    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--cell-type-twas", action="store_true",
                      help="run the cell-type TWAS pipeline (default when no step selector is given)")
    mode.add_argument("--all", dest="cell_type_twas", action="store_true", help=argparse.SUPPRESS)
    mode.add_argument("--steps", nargs="+", metavar="STEP",
                      help="run only these steps in pipeline order: " + ",".join(ALL_STEPS))

    p.add_argument("--chromosomes", default=",".join(str(c) for c in range(1, 23)),
                   help="comma-separated chromosome list (default 1-22)")
    p.add_argument("--cell-types", required=True,
                   help="space/comma-separated cell types matching the reference and fractions")

    data = p.add_argument_group("data inputs")
    data.add_argument("--bulk-tpm", default="", help="bulk TPM expression matrix (genes x samples, tab-separated)")
    data.add_argument("--enigma-ref", default="", help="single-cell reference matrix (TPM, aggregated by cell type)")
    data.add_argument("--cibersortx-frax", default="", help="CIBERSORTx cell-fraction output passed to ENIGMA")
    data.add_argument("--covariates", default="", help="covariates txt (required by peerCov)")
    data.add_argument("--snp-annot-pattern", default="", help="per-chromosome SNP annotation, {chr} placeholder")
    data.add_argument("--genotype-pattern", default="", help="per-chromosome genotype dosage, {chr} placeholder")
    data.add_argument("--gene-annot", default="", help="gene annotation file (e.g. gene_annot.parsed.txt)")
    data.add_argument("--snp-annot-file", default="", help="single-file genotype panel SNP annotation (used by alignGWAS)")
    data.add_argument("--gwas-dir", default="", help="directory holding GWAS files (used by legacy --gwas-profile)")
    data.add_argument("--g1000-vcf-pattern", default="", help="1000G Phase 3 VCF pattern with {chr} placeholder")
    data.add_argument("--g1000-samples", default="", help="1000G sample list matching the GWAS ancestry")
    data.add_argument("--gene-pos", default="", help="gene position annotation CSV (gene_id/gene_name columns)")
    data.add_argument("--ld-dosage-dir", default="", help="pre-extracted 1000G dosage directory (skips extractLD extraction)")
    data.add_argument("--gwas-summary", default="",
                      help="pre-QC'd GWAS summary for runTWAS (skips alignGWAS/metaGWAS/filterGWAS)")

    gen = p.add_argument_group("input contracts and generalization")
    gen.add_argument("--sample-id-policy", choices=["exact", "legacy-rnames"], default=None,
                     help="sample ID matching between expression and covariate/genotype files (default exact)")
    gen.add_argument("--cov-layout", choices=["auto", "samples", "covariates"], default="auto",
                     help="peerCov covariates file layout (default auto)")
    gen.add_argument("--gene-annot-format", choices=["standard", "legacy-name-is-id"], default=None,
                     help="gene annotation schema (default standard)")
    gen.add_argument("--genome-build", choices=["GRCh37", "GRCh38"], default=None,
                     help="genome build declared for this cohort's coordinates")
    gen.add_argument("--training-population", default="unspecified",
                     help="optional population label stored with the models (default unspecified)")
    gen.add_argument("--maf", type=float, default=0.01,
                     help="MAF threshold for trainModel genotypes in [0, 0.5] (default 0.01)")
    gen.add_argument("--maf-policy", choices=["filter", "verify"], default="filter",
                     help="filter: drop variants < --maf; verify: fail if genotype violates --maf (default filter)")
    gen.add_argument("--gwas-population", default=None,
                     help="GWAS ancestry label (e.g. EUR, EAS, EUR_EAS); defaults --ld-panel-id when omitted")
    gen.add_argument("--ld-panel-id", default=None,
                     help="LD panel id matched to GWAS ancestry (defaults to --gwas-population, else EUR_EAS)")
    gen.add_argument("--twas-model-set", choices=["unfiltered", "filtered-signif"], default="filtered-signif",
                     help="weight DB for runTWAS (default filtered-signif; buildDB always retains unfiltered DB for buildLD)")
    gen.add_argument("--cores", type=int, default=None, help="parallel cores/workers across pipeline steps (default 4)")
    gen.add_argument("--seed", type=int, default=None, help="explicit trainModel seed; forces --training-cores 1")
    gen.add_argument("--training-cores", type=int, default=None, help="inner trainModel cores per (cell_type, chr) job")

    gwas = p.add_argument_group("GWAS mode (required when running alignGWAS/metaGWAS/filterGWAS)")
    gwas_mode = gwas.add_mutually_exclusive_group()
    gwas_mode.add_argument("--gwas-profile", choices=list(LEGACY_GWAS_PROFILES), default=None,
                           help="historical SLE reproduction profile (retained for backward compatibility)")
    gwas_mode.add_argument("--gwas-manifest", default=None, metavar="FILE",
                           help="TSV manifest of GWAS studies (see README for column contract)")
    gwas.add_argument("--gwas-filter-policy", choices=["none", "legacy"], default=None,
                      help="post-meta filtering (default: none for --gwas-manifest, legacy for --gwas-profile)")
    gwas.add_argument("--skip-gwas-meta", action="store_true",
                      help="skip optional metaGWAS/filterGWAS for a single pre-QC'd study")

    env = p.add_argument_group("environment and tools")
    env.add_argument("--env-base", default="", help="portable conda-pack environment root")
    env.add_argument("--r-env", default="R4.4.3", help="conda env with R 4.4.3 + ENIGMA")
    env.add_argument("--peer-env", default="r3.6.3", help="conda env with R 3.6.3 + PEER")
    env.add_argument("--py-env", default="imlatools", help="conda env with MetaXcan runtime")
    env.add_argument("--metaxcan-dir", default=os.path.join(REPO, "tools", "MetaXcan"))
    env.add_argument("--metal-bin", default=os.path.join(REPO, "tools", "METAL", "metal"))
    env.add_argument("--bcftools", default="bcftools", help="bcftools executable")
    env.add_argument("--bcf-plugins", default="", help="bcftools plugins directory (empty = auto)")
    env.add_argument("--gwas-id", default=None, help="trait label used in output names and reports")

    step_flags = p.add_argument_group("step shortcuts (alternative to --cell-type-twas / --steps)")
    for step in ALL_STEPS:
        step_flags.add_argument("--" + step, dest="flag_" + step, action="store_true",
                                help="run " + step + ": " + STEP_DESCRIPTIONS[step])
    return p


def _default_twas_steps(args):
    if args.gwas_summary:
        return [s for s in ALL_STEPS if s not in GWAS_STEPS]
    if args.skip_gwas_meta:
        return [s for s in ALL_STEPS if s not in ("metaGWAS", "filterGWAS")]
    return ALL_STEPS.copy()


def _validate_gwas_inputs(parser, args):
    """Validate direct --gwas-summary and --skip-gwas-meta flag combinations."""
    if args.gwas_summary:
        if args.gwas_profile or args.gwas_manifest:
            parser.error("--gwas-summary cannot be combined with --gwas-profile or --gwas-manifest")
        if args.skip_gwas_meta:
            parser.error("--skip-gwas-meta cannot be combined with --gwas-summary")
        if set(args.steps) & GWAS_STEPS:
            parser.error("--gwas-summary cannot be combined with alignGWAS, metaGWAS, or filterGWAS")
        if not os.path.isfile(args.gwas_summary) or os.path.getsize(args.gwas_summary) == 0:
            parser.error("--gwas-summary file not found or empty: " + args.gwas_summary)
        opener = gzip.open if args.gwas_summary.endswith(".gz") else open
        with opener(args.gwas_summary, "rt") as gf:
            hdr = gf.readline().rstrip("\r\n").split()
        if not (GWAS_SUMMARY_REQUIRED_COLS <= set(hdr)):
            parser.error("--gwas-summary header must contain columns: " + ", ".join(sorted(GWAS_SUMMARY_REQUIRED_COLS)))

    if args.skip_gwas_meta:
        if args.gwas_profile:
            parser.error("--skip-gwas-meta applies only to single-study --gwas-manifest (not --gwas-profile)")
        if set(args.steps) & {"metaGWAS", "filterGWAS"}:
            parser.error("--skip-gwas-meta cannot be combined with metaGWAS or filterGWAS")


def _recover_gwas_mode(parser, args, saved, run_meta, needs_mode):
    """Validate and recover prior GWAS configuration from gwas_run.json."""
    if not isinstance(saved, dict):
        parser.error("Invalid GWAS mode record: " + run_meta)
    saved_id = saved.get("gwas_id")
    if not isinstance(saved_id, str) or not re.fullmatch(SAFE_NAME, saved_id):
        parser.error("Invalid GWAS identity record: " + run_meta)
    if args.gwas_id is not None and args.gwas_id != saved_id:
        parser.error("--gwas-id conflicts with recorded GWAS inputs; use a new output directory")

    profile, manifest = args.gwas_profile, args.gwas_manifest
    if not needs_mode:
        print(f"[gwas] recovered gwas_id={saved_id!r} from {run_meta}")
        args.gwas_id = saved_id
        return profile, manifest

    req_profile, req_manifest = profile, manifest
    profile, manifest = saved.get("profile"), saved.get("manifest")
    if (
        bool(profile) == bool(manifest)
        or (profile and profile not in LEGACY_GWAS_PROFILES)
        or (manifest and (not isinstance(manifest, str) or not os.path.isabs(manifest)))
    ):
        parser.error("Invalid GWAS mode record: " + run_meta)
    if (req_profile and req_profile != profile) or (
        req_manifest and (not manifest or os.path.abspath(req_manifest) != manifest)
    ):
        parser.error("Explicit GWAS mode conflicts with existing provenance; use a new output directory")
    if manifest and saved.get("manifest_sha256"):
        if os.path.isfile(manifest):
            if sha256_file(manifest) != saved["manifest_sha256"]:
                parser.error("GWAS manifest content changed; use a new output directory")
        elif req_manifest:
            parser.error("Explicit GWAS manifest no longer exists")

    saved_build = saved.get("genome_build", "GRCh37")
    saved_policy = saved.get("filter_policy")
    saved_pop = saved.get("gwas_population") or ""
    if args.gwas_population is not None and saved_pop and args.gwas_population != saved_pop:
        parser.error("--gwas-population conflicts with recorded GWAS inputs; use a new output directory")
    if args.gwas_population is None and saved_pop:
        args.gwas_population = saved_pop
    if args.skip_gwas_meta and not bool(saved.get("skip_meta", False)):
        parser.error("--skip-gwas-meta conflicts with recorded meta-analysis mode; use a new output directory")
    if (
        saved_build not in ("GRCh37", "GRCh38")
        or saved_policy not in ("none", "legacy")
        or (profile and saved_policy != "legacy")
    ):
        parser.error("Invalid GWAS identity/build/filter record: " + run_meta)
    for flag, old, new in (("--genome-build", saved_build, args.genome_build),
                           ("--gwas-filter-policy", saved_policy, args.gwas_filter_policy)):
        if new is not None and new != old:
            parser.error(f"{flag} conflicts with recorded GWAS inputs; use a new output directory")

    args.genome_build = saved_build
    args.gwas_filter_policy = saved_policy
    args.gwas_mode_recovered = True
    args.gwas_id = saved_id
    print("[gwas] recovered recorded mode from", run_meta)
    return profile, manifest


def resolve_twas_gwas_mode(parser, args):
    """Resolve GWAS input mode, recover provenance if cached, and set output paths."""
    _validate_gwas_inputs(parser, args)

    gwas_dir = os.path.join(os.path.abspath(args.output_dir), "M04_gwas")
    run_meta = os.path.join(gwas_dir, "gwas_run.json")
    needs_mode = bool(set(args.steps) & GWAS_STEPS) or ("runTWAS" in args.steps and not args.gwas_summary)
    needs_summary_id = "summarizeTWAS" in args.steps
    args.gwas_mode_recovered = False
    profile, manifest = args.gwas_profile, args.gwas_manifest
    saved = None

    if needs_mode or (needs_summary_id and not args.gwas_summary):
        saved = read_metadata(run_meta)
        if saved is not None:
            profile, manifest = _recover_gwas_mode(parser, args, saved, run_meta, needs_mode)

    if needs_mode and not profile and not manifest:
        parser.error("GWAS steps require --gwas-manifest or --gwas-profile, or a recorded mode in " + run_meta)

    if args.gwas_id is None:
        if profile:
            args.gwas_id = "SLE"
        elif manifest:
            parser.error("--gwas-manifest requires an explicit --gwas-id")
        elif needs_summary_id:
            parser.error("summarizeTWAS requires --gwas-id when no recorded GWAS provenance exists in " + run_meta)
        else:
            args.gwas_id = ""

    args.genome_build = args.genome_build or "GRCh37"
    if profile:
        if args.genome_build != "GRCh37":
            parser.error("Historical SLE profiles require --genome-build GRCh37")
        args.gwas_prefix, args.skip_meta_active, filter_policy = LEGACY_GWAS_PREFIX, False, "legacy"
    elif manifest:
        is_single = (
            args.skip_gwas_meta
            or (saved is not None and bool(saved.get("skip_meta", False)))
            or ("alignGWAS" in args.steps and "metaGWAS" not in args.steps and "filterGWAS" not in args.steps)
        )
        if is_single and (not args.gwas_mode_recovered or os.path.isfile(manifest)):
            with open(manifest, "r", newline="") as mf:
                rows = list(csv.DictReader(mf, delimiter="\t"))
            if len(rows) != 1:
                parser.error(f"skipping metaGWAS/filterGWAS requires a single-study --gwas-manifest (found {len(rows)} studies)")
        args.gwas_prefix = args.gwas_id + ("_single" if is_single else "_meta")
        args.skip_meta_active = is_single
        filter_policy = "none"
    else:
        args.gwas_prefix, args.skip_meta_active, filter_policy = "", False, ""

    if args.gwas_filter_policy is not None:
        if profile and args.gwas_filter_policy != "legacy":
            parser.error("the legacy profile fixes --gwas-filter-policy legacy")
        filter_policy = args.gwas_filter_policy
    if needs_mode and saved and saved.get("prefix") != args.gwas_prefix:
        parser.error("Recorded GWAS prefix disagrees with analysis identity")

    args.gwas_filter_policy = filter_policy
    args.metal_tbl = os.path.join(gwas_dir, args.gwas_prefix + "1.tbl") if args.gwas_prefix else ""
    args.gwas_filt = (
        os.path.abspath(args.gwas_summary)
        if args.gwas_summary
        else (os.path.join(gwas_dir, args.gwas_prefix + "_filt.tsv") if args.gwas_prefix else "")
    )
    args.ld_panel_id = args.ld_panel_id or args.gwas_population or "EUR_EAS"
    args.gwas_profile_active = profile
    args.gwas_manifest_active = manifest


def parse_args(argv=None):
    """Parse and validate arguments for the cell-type TWAS pipeline."""
    p = build_twas_parser()
    args = p.parse_args(argv)

    args.steps, args.step_aliases = resolve_step_selection(
        p, args, ALL_STEPS, STEP_ALIASES, _default_twas_steps,
        empty_msg="--steps cannot contain empty steps",
        unknown_prefix="unknown steps",
        has_full_flag=True,
    )

    check_finite_bounds(p, args.alpha, "--alpha", gt=0)
    check_finite_bounds(p, args.min_pct, "--min-pct", ge=0, le=1)
    check_finite_bounds(p, args.tpm_min, "--tpm-min", ge=0)
    check_finite_bounds(p, args.min_samples, "--min-samples", ge=0)
    check_finite_bounds(p, args.min_level_count, "--min-level-count", ge=0)
    check_finite_bounds(p, args.maf, "--maf", ge=0, le=0.5)
    if args.seed is not None and not 0 <= args.seed <= 2147483647:
        p.error("--seed must be an integer in 0..2147483647")
    if args.seed is not None and args.training_cores not in (None, 1):
        p.error("--seed requires serial inner training; set --training-cores 1 or omit it")

    args.cores_resolved = resolve_cores(p, args.cores, args.training_cores, "--training-cores")
    args.chromosomes = parse_chromosomes(p, args.chromosomes)
    for flag, value in (("--gwas-id", args.gwas_id), ("--gwas-population", args.gwas_population), ("--ld-panel-id", args.ld_panel_id)):
        validate_safe_label(p, flag, value)
    args.cell_types = parse_cell_types(p, args.cell_types, required=True)

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
        p.error("--twas-model-set (unfiltered or filtered-signif) is required for buildDB and runTWAS")

    required_inputs = {
        "inferExp": ("bulk_tpm", "enigma_ref", "cibersortx_frax"),
        "peerCov": ("covariates",),
        "trainModel": ("snp_annot_pattern", "genotype_pattern", "gene_annot"),
        "alignGWAS": ("snp_annot_file",) if args.gwas_manifest else ("snp_annot_file", "gwas_dir"),
        "extractLD": ("snp_annot_pattern",) if args.ld_dosage_dir else ("snp_annot_pattern", "g1000_vcf_pattern", "g1000_samples"),
        "summarizeTWAS": ("gene_pos",),
    }
    missing = sorted({
        "--" + name.replace("_", "-")
        for step in args.steps
        for name in required_inputs.get(step, ())
        if not getattr(args, name)
    })
    if missing:
        p.error("selected steps require: " + ", ".join(missing))

    resolve_twas_gwas_mode(p, args)
    return args


# ---------------------------------------------------------------------------
# LINCS therapeutic prioritization parser and runner
# ---------------------------------------------------------------------------

def build_lincs_parser():
    """Construct the argparse parser for `medicell.py lincs`."""
    p = argparse.ArgumentParser(
        prog="medicell.py lincs",
        description="MediCell LINCS therapeutic prioritization: perturbational enrichment and cell-type scoring",
    )
    p.add_argument("--output-dir", required=True,
                   help="output directory (may reuse TWAS output-dir containing M01_deconvolution/M08_summary)")
    p.add_argument("--steps", nargs="+", metavar="STEP", default=None,
                   help="LINCS steps in order: " + ",".join(LINCS_STEPS))
    for lstep in LINCS_STEPS:
        p.add_argument("--" + lstep, dest="flag_" + lstep, action="store_true", help="run LINCS step " + lstep)
    p.add_argument("--trait", default="", help="trait label for output tags (default: from M04 gwas_run.json or 'trait')")
    p.add_argument("--tissue", "--label", dest="tissue", default="cohort", help="cohort label in output names (default cohort)")
    p.add_argument("--direction", choices=["pos", "neg", "both"], default="pos", help="TWAS z-score direction (default pos)")
    p.add_argument("--cell-types", default="", help="optional space/comma-separated cell types to retain")

    twas = p.add_argument_group("TWAS and expression inputs")
    twas.add_argument("--twas-file", default="", help="TWAS summary table (.csv/.tsv/.xlsx)")
    twas.add_argument("--twas-sheet", default="", help="Excel sheet name when --twas-file is .xlsx")
    twas.add_argument("--twas-gene-col", default="", help="gene symbol column name in --twas-file")
    twas.add_argument("--twas-zscore-col", default="", help="z-score column name in --twas-file")
    twas.add_argument("--twas-padj-col", default="", help="adjusted p-value column name in --twas-file")
    twas.add_argument("--twas-celltype-col", default="", help="cell type column name in --twas-file")
    twas.add_argument("--twas-padj-cutoff", type=float, default=0.05, help="TWAS adjusted p-value cutoff (default 0.05)")
    twas.add_argument("--twas-zscore-cutoff", type=float, default=0.25, help="minimum |z-score| cutoff (default 0.25)")
    twas.add_argument("--exp-matrix", default="", help="pre-built gene x cell-type expression CSV")
    twas.add_argument("--enigma-rdata", default="", help="inferExp ENIGMA RData path")
    twas.add_argument("--alpha", type=float, default=0.1, help="ENIGMA alpha for default RData path (default 0.1)")
    twas.add_argument("--gene-pos", default="", help="gene position CSV (gene_id, gene_name columns)")

    lincs = p.add_argument_group("LINCS enrichment and scoring")
    lincs.add_argument("--lincs-gmt", default=os.path.join(REPO, "data", "ref", "lincs", "LINCS_PRL_top250_human.gmt"))
    lincs.add_argument("--lincs-pert-phase1", default=os.path.join(REPO, "data", "ref", "lincs", "GSE92742_Broad_LINCS_pert_info.txt.gz"))
    lincs.add_argument("--lincs-pert-phase2", default=os.path.join(REPO, "data", "ref", "lincs", "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"))
    lincs.add_argument("--lincs-input-csv", default="", help="pre-split directional TWAS CSV for standalone enrichLincs")
    lincs.add_argument("--lincs-sig-csv", default="", help="pre-computed LINCS enrichment CSV for scoreLincs")
    lincs.add_argument("--lincs-min-size", type=int, default=10, help="fgsea minSize (default 10)")
    lincs.add_argument("--lincs-max-size", type=int, default=500, help="fgsea maxSize (default 500)")
    lincs.add_argument("--lincs-sig-col", choices=["pval", "padj"], default="pval", help="significance filter column (default pval)")
    lincs.add_argument("--lincs-sig-cutoff", type=float, default=0.05, help="significance cutoff (default 0.05)")
    lincs.add_argument("--lincs-nes-filter", choices=["neg", "pos", "all"], default="neg", help="NES filter (default neg)")
    lincs.add_argument("--cores", type=int, default=None, help="parallel cores across LINCS steps (default 4)")
    lincs.add_argument("--lincs-cores", type=int, default=None, help="fgsea thread override (defaults to --cores)")
    lincs.add_argument("--remove-ribo", action="store_true", help="exclude ribosomal/mitochondrial genes before fgsea")
    lincs.add_argument("--top-n", type=int, default=50, help="top perturbations retained per run (default 50)")
    lincs.add_argument("--lincs-rank-by", choices=["pval", "padj"], default="pval", help="ranking metric in scoreLincs (default pval)")
    lincs.add_argument("--exclude-pert-pattern", default=r"^BRD-|trt_sh\.css", help="regex to exclude perturbations ('' disables)")

    env = p.add_argument_group("environment")
    env.add_argument("--env-base", default="", help="portable conda-pack environment root")
    env.add_argument("--lincs-env", default="", help="conda env with fgsea/dplyr/tidyr/readxl (default: auto-detect)")
    return p


def parse_lincs_args(argv):
    """Parse and validate arguments for `medicell.py lincs`."""
    p = build_lincs_parser()
    args = p.parse_args(argv)

    validate_safe_label(p, "--trait", args.trait)
    validate_safe_label(p, "--tissue/--label", args.tissue)
    args.steps, _ = resolve_step_selection(
        p, args, LINCS_STEPS, LINCS_STEP_ALIASES,
        lambda a: ["prepLincs", "scoreLincs"] if a.lincs_sig_csv else LINCS_STEPS.copy(),
        empty_msg="--steps cannot contain empty values",
        unknown_prefix="unknown LINCS steps",
        has_full_flag=False,
    )
    args.cell_types = parse_cell_types(p, args.cell_types, required=False)

    if "scoreLincs" in args.steps and len(args.cell_types) == 1:
        p.error("LINCS scoring (scoreLincs) requires at least 2 cell types; single-cell-type input yields undefined cross-cell-type expression z-scores")
    check_finite_bounds(p, args.top_n, "--top-n", ge=1)
    if args.lincs_min_size < 1 or args.lincs_max_size < args.lincs_min_size:
        p.error("--lincs-min-size and --lincs-max-size must satisfy 1 <= min <= max")

    args.cores_resolved = resolve_cores(p, args.cores, args.lincs_cores, "--lincs-cores")
    args.lincs_cores_resolved = args.lincs_cores if args.lincs_cores is not None else args.cores_resolved
    if not (0 < args.twas_padj_cutoff <= 1) or not (0 < args.lincs_sig_cutoff <= 1):
        p.error("p-value cutoffs must be in (0, 1]")
    check_finite_bounds(p, args.twas_zscore_cutoff, "--twas-zscore-cutoff", ge=0)
    return args


def _resolve_lincs_inputs(args, out):
    """Resolve trait label, TWAS table, and ENIGMA expression matrix paths."""
    trait = args.trait
    if not trait:
        gwas_meta = read_metadata(os.path.join(out, "M04_gwas", "gwas_run.json"))
        saved_id = gwas_meta.get("gwas_id") if isinstance(gwas_meta, dict) else None
        trait = saved_id if isinstance(saved_id, str) and saved_id else "trait"

    twas_file = abspath_opt(args.twas_file)
    if not twas_file:
        for fname in ("allcelltype_result.csv", "FDRsig_allcelltype_result.csv"):
            candidate = os.path.join(out, "M08_summary", fname)
            if os.path.isfile(candidate):
                twas_file = candidate
                break

    enigma_rdata = abspath_opt(args.enigma_rdata)
    if not enigma_rdata and not args.exp_matrix:
        candidate = os.path.join(out, "M01_deconvolution", f"enigma.cse.trace_alpha_{args.alpha}.RData")
        if os.path.isfile(candidate):
            enigma_rdata = candidate
    return trait, twas_file, enigma_rdata


def build_lincs_env(args, dirs, trait, twas_file, enigma_rdata, direction, step_cores):
    """Build MC_* environment variables for one LINCS direction run."""
    twas_tag = f"{trait}_{args.tissue}_{direction}"
    l01_dir_csv = os.path.join(dirs["l01"], f"{twas_tag}.csv")
    l01_twas_tsv = os.path.join(dirs["l01"], f"{trait}_{args.tissue}_twas_all.tsv")
    l01_exp_csv = os.path.join(dirs["l01"], f"{args.tissue}_exp_matrix.csv")
    l02_sig_csv = os.path.join(dirs["l02"], f"LINCS_sig_{twas_tag}.csv")

    l02_in_csv = abspath_opt(args.lincs_input_csv) if args.lincs_input_csv else l01_dir_csv
    l03_sig_csv = abspath_opt(args.lincs_sig_csv) if args.lincs_sig_csv else l02_sig_csv
    l03_twas_tsv = l01_twas_tsv if ("prepLincs" in args.steps or os.path.isfile(l01_twas_tsv)) else twas_file
    l03_exp_csv = l01_exp_csv if ("prepLincs" in args.steps or os.path.isfile(l01_exp_csv)) else abspath_opt(args.exp_matrix)

    if "scoreLincs" in args.steps and (not l03_twas_tsv or not l03_exp_csv):
        sys.exit("scoreLincs requires prepLincs outputs or explicit --twas-file and --exp-matrix")

    use_l01_twas = l03_twas_tsv == l01_twas_tsv
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
        "MC_LINCS_L03_GENE_COL": "" if use_l01_twas else args.twas_gene_col,
        "MC_LINCS_L03_CT_COL": "" if use_l01_twas else args.twas_celltype_col,
        "MC_LINCS_L03_Z_COL": "" if use_l01_twas else args.twas_zscore_col,
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
        "MC_LINCS_EXCLUDE_PATTERN": args.exclude_pert_pattern or "NONE",
        **blas_env(1),
    })
    return env


def run_lincs_mode(argv):
    """Execute selected LINCS therapeutic prioritization steps."""
    args = parse_lincs_args(argv)
    env_base = abspath_opt(args.env_base)
    lincs_env = resolve_lincs_env(env_base, args.lincs_env)
    conda_run = make_conda_runner(env_base)

    out = os.path.abspath(args.output_dir)
    dirs = ensure_dirs({
        "l01": os.path.join(out, "L01_inputs"),
        "l02": os.path.join(out, "L02_enrichment"),
        "l03": os.path.join(out, "L03_drug_scores"),
    })
    trait, twas_file, enigma_rdata = _resolve_lincs_inputs(args, out)

    if "prepLincs" in args.steps:
        if not twas_file:
            sys.exit("prepLincs requires --twas-file (or M08_summary/allcelltype_result.csv in --output-dir)")
        if not args.exp_matrix and not enigma_rdata:
            sys.exit("prepLincs requires --exp-matrix or --enigma-rdata (or M01_deconvolution RData in --output-dir)")
    if "enrichLincs" in args.steps:
        if not os.path.isfile(args.lincs_gmt):
            sys.exit(f"enrichLincs requires existing file for --lincs-gmt: {args.lincs_gmt}")
        for flag, fpath in (("--lincs-pert-phase1", args.lincs_pert_phase1),
                            ("--lincs-pert-phase2", args.lincs_pert_phase2)):
            if fpath and not os.path.isfile(fpath):
                sys.exit(f"enrichLincs annotation file not found for {flag}: {fpath}")

    directions = ["pos", "neg"] if args.direction == "both" else [args.direction]
    if args.direction == "both" and (args.lincs_sig_csv or args.lincs_input_csv):
        sys.exit("--lincs-sig-csv / --lincs-input-csv cannot be combined with --direction both")

    total_cores = args.cores_resolved
    dir_workers = min(len(directions), total_cores)
    l02_inner_cores = max(1, args.lincs_cores_resolved // dir_workers)

    print("[lincs steps]", " -> ".join(args.steps))
    print("[lincs trait/label]", trait, "/", args.tissue, "directions:", ",".join(directions), "cores:", total_cores)
    print("[lincs env]", lincs_env)

    step_workers_cores = {
        "prepLincs": (1, total_cores),
        "enrichLincs": (dir_workers, l02_inner_cores),
        "scoreLincs": (dir_workers, 1),
    }
    for step in args.steps:
        workers, step_cores = step_workers_cores[step]
        cmd = script_cmd(os.path.join(SCRIPTS, LINCS_SCRIPT_FILES[step]))
        run_parallel(
            workers,
            directions,
            lambda d, c=cmd, sc=step_cores: conda_run(
                lincs_env, c, build_lincs_env(args, dirs, trait, twas_file, enigma_rdata, d, sc)
            ),
        )

    print("MediCell LINCS run complete. Results in:", out)


# ---------------------------------------------------------------------------
# Cell-type TWAS execution and CLI dispatch
# ---------------------------------------------------------------------------

def resolve_m03_parallelism(total_cores, n_jobs, seed, training_cores):
    """Allocate outer (cell_type, chr) workers and inner OpenMP/doParallel cores."""
    if seed is not None:
        return min(n_jobs, total_cores), 1
    if training_cores is not None:
        return max(1, min(n_jobs, total_cores // training_cores)), training_cores
    workers = min(n_jobs, total_cores)
    return workers, max(1, total_cores // workers)


def build_twas_base_env(args, dirs):
    """Build shared MC_* environment variables passed to all TWAS steps."""
    chrs = args.chromosomes
    return {
        "MC_ALPHA": str(args.alpha),
        "MC_TPM_MIN": str(args.tpm_min),
        "MC_MIN_PCT": str(args.min_pct),
        "MC_MIN_SAMPLES": str(args.min_samples),
        "MC_BULK_TPM": abspath_opt(args.bulk_tpm),
        "MC_ENIGMA_REF": abspath_opt(args.enigma_ref),
        "MC_CBS_FRAX": abspath_opt(args.cibersortx_frax),
        "MC_COVARIATES": abspath_opt(args.covariates),
        "MC_CAT_COVS": args.cat_covs,
        "MC_MIN_LEVEL_COUNT": str(args.min_level_count),
        "MC_PEER_NK": "" if args.peer_nk_resolved is None else str(args.peer_nk_resolved),
        "MC_PEER_KEEP": "" if args.peer_keep_resolved is None else str(args.peer_keep_resolved),
        "MC_SNP_ANNOT_PATTERN": abspath_opt(args.snp_annot_pattern),
        "MC_GENOTYPE_PATTERN": abspath_opt(args.genotype_pattern),
        "MC_GENE_ANNOT": abspath_opt(args.gene_annot),
        "MC_GENE_ANNOT_FORMAT": args.gene_annot_format,
        "MC_SNP_ANNOT": abspath_opt(args.snp_annot_file),
        "MC_GWAS_IN_DIR": abspath_opt(args.gwas_dir),
        "MC_GWAS_PROFILE": args.gwas_profile_active or "",
        "MC_GWAS_MANIFEST": abspath_opt(args.gwas_manifest_active),
        "MC_GWAS_ID": args.gwas_id,
        "MC_GWAS_PREFIX": args.gwas_prefix,
        "MC_GWAS_OUT": dirs["gwas"],
        "MC_GWAS_FILTER_POLICY": args.gwas_filter_policy,
        "MC_METAL_TBL": args.metal_tbl,
        "MC_GWAS_FILT": args.gwas_filt,
        "MC_GWAS_SINGLE_OUT": args.gwas_filt if args.skip_meta_active else "",
        "MC_GWAS_POPULATION": args.gwas_population or "",
        "MC_G1000_VCF_PATTERN": abspath_opt(args.g1000_vcf_pattern),
        "MC_G1000_SAMPLES": abspath_opt(args.g1000_samples),
        "MC_GENE_POS": abspath_opt(args.gene_pos),
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
        "MC_TRAINING_CORES": str(1 if args.seed is not None else (args.training_cores or args.cores_resolved)),
        "MC_METAXCAN_DIR": os.path.abspath(args.metaxcan_dir),
        "MC_METAL_BIN": os.path.abspath(args.metal_bin),
        "MC_BCFTOOLS": args.bcftools,
        "MC_BCF_PLUGINS": args.bcf_plugins,
        "MC_ELNET_SCRIPT": os.path.join(SCRIPTS, "helpers", "gtex_v7_nested_cv_elnet.R"),
        "MC_CELLTYPES": " ".join(args.cell_types),
        "MC_CHROMS": " ".join(str(c) for c in chrs),
        "MC_CHROMS_CSV": ",".join(str(c) for c in chrs),
        "MC_LD_SKIP_EXTRACT": "1" if args.ld_dosage_dir else "0",
    }


def record_twas_gwas_mode(args, gwas_dir):
    """Persist GWAS provenance to M04_gwas/gwas_run.json on first execution."""
    if args.gwas_mode_recovered:
        return
    digest = (
        sha256_file(args.gwas_manifest_active)
        if args.gwas_manifest_active and os.path.isfile(args.gwas_manifest_active)
        else ""
    )
    write_metadata(os.path.join(gwas_dir, "gwas_run.json"), {
        "manifest_sha256": digest,
        "profile": args.gwas_profile_active or "",
        "manifest": abspath_opt(args.gwas_manifest_active),
        "gwas_id": args.gwas_id,
        "genome_build": args.genome_build,
        "prefix": args.gwas_prefix,
        "metal_tbl": args.metal_tbl,
        "filtered": args.gwas_filt,
        "filter_policy": args.gwas_filter_policy,
        "gwas_population": args.gwas_population or "",
        "skip_meta": bool(args.skip_meta_active),
    })


def run_twas_mode(argv):
    """Execute selected cell-type TWAS steps in dependency order."""
    args = parse_args(argv)
    env_base = abspath_opt(args.env_base)
    if env_base and args.metal_bin == os.path.join(REPO, "tools", "METAL", "metal"):
        args.metal_bin = os.path.join(env_base, "r3.6.3", "bin", "metal")
    conda_run = make_conda_runner(env_base)

    if args.step_aliases:
        print("[step aliases]", ", ".join(f"{old}->{canon}" for old, canon in args.step_aliases))
    print("[steps]", " -> ".join(args.steps))
    if args.gwas_profile_active:
        print(f"[gwas] WARNING: legacy profile {args.gwas_profile_active!r} reproduces historical SLE behavior")

    out = os.path.abspath(args.output_dir)
    dirs = ensure_dirs({
        "m01": os.path.join(out, "M01_deconvolution"),
        "exp": os.path.join(out, "M02_expression"),
        "cov": os.path.join(out, "M02_covariates"),
        "m03": os.path.join(out, "M03_models"),
        "db": os.path.join(out, "M06_databases"),
        "gwas": os.path.join(out, "M04_gwas"),
        "ld_dosage": args.ld_dosage_dir or os.path.join(out, "M05_ld", "dosage"),
        "ld_cov": os.path.join(out, "M05_ld", "cov_" + args.ld_panel_id),
        "twas": os.path.join(out, "M07_twas"),
        "summary": os.path.join(out, "M08_summary"),
    })

    base = build_twas_base_env(args, dirs)

    def mc(**overrides):
        env = os.environ.copy()
        env.update(base)
        env.update(overrides)
        return env

    def twas_cmd(step, *extra_args):
        return script_cmd(os.path.join(SCRIPTS, SCRIPT_FILES[step]), *extra_args)

    total_cores = args.cores_resolved
    single_blas, multi_blas = blas_env(1), blas_env(total_cores)
    enigma_rdata = os.path.join(dirs["m01"], f"enigma.cse.trace_alpha_{args.alpha}.RData")
    save_gwas = lambda: record_twas_gwas_mode(args, dirs["gwas"])

    step_plan = {
        "inferExp": (args.r_env, {"MC_OUT_DIR": dirs["m01"], **multi_blas}, None),
        "normExp": (args.peer_env, {"MC_ENIGMA_RDATA": enigma_rdata, "MC_EXP_DIR": dirs["exp"], **single_blas}, None),
        "peerCov": (args.peer_env, {"MC_EXP_DIR": dirs["exp"], "MC_COV_DIR": dirs["cov"], **single_blas}, None),
        "buildDB": (args.r_env, {
            "MC_M03_SUMMARY_DIR": os.path.join(dirs["m03"], "summary"),
            "MC_M03_WEIGHTS_DIR": os.path.join(dirs["m03"], "weights"),
            "MC_DB_DIR": dirs["db"],
            "MC_CHROMS": base["MC_CHROMS_CSV"],
            **single_blas,
        }, None),
        "alignGWAS": (args.r_env, {}, save_gwas),
        "metaGWAS": (args.r_env, {"MC_GWAS_OUT": dirs["gwas"], "MC_METAL_TBL": args.metal_tbl}, save_gwas),
        "filterGWAS": (args.r_env, {
            "MC_METAL_TBL": args.metal_tbl, "MC_GWAS_OUT": dirs["gwas"], "MC_GWAS_FILT": args.gwas_filt,
        }, save_gwas),
        "extractLD": (args.py_env, {"MC_LD_DOSAGE_DIR": dirs["ld_dosage"], "MC_LD_COV_DIR": dirs["ld_cov"], **single_blas}, None),
        "buildLD": (args.py_env, {"MC_DB_DIR": dirs["db"], "MC_LD_COV_DIR": dirs["ld_cov"], **single_blas}, None),
        "runTWAS": (args.py_env, {
            "MC_DB_DIR": dirs["db"], "MC_LD_COV_DIR": dirs["ld_cov"],
            "MC_GWAS_FILT": args.gwas_filt, "MC_TWAS_OUT": dirs["twas"], **single_blas,
        }, None),
        "summarizeTWAS": (args.r_env, {
            "MC_TWAS_DIR": dirs["twas"], "MC_SUMMARY_OUT": dirs["summary"],
            "MC_GWAS_ID": args.gwas_id, **single_blas,
        }, None),
    }

    for step in args.steps:
        if step == "trainModel":
            jobs = [(ct, chr_) for ct in args.cell_types for chr_ in args.chromosomes]
            workers, inner_cores = resolve_m03_parallelism(total_cores, len(jobs), args.seed, args.training_cores)
            run_parallel(
                workers,
                jobs,
                lambda job: conda_run(
                    args.r_env,
                    twas_cmd("trainModel", job[0], job[1]),
                    mc(
                        MC_COV_FILE=os.path.join(dirs["cov"], f"{job[0]}_covariates.txt"),
                        MC_EXP_FILE=os.path.join(dirs["exp"], f"{job[0]}_transformed_expression.txt"),
                        MC_M03_OUT=dirs["m03"],
                        MC_TRAINING_CORES=str(inner_cores),
                        **single_blas,
                    ),
                ),
            )
            continue

        env_name, overrides, post_hook = step_plan[step]
        conda_run(env_name, twas_cmd(step), mc(**overrides))
        if post_hook is not None:
            post_hook()

    print("MediCell run complete. Results in:", out)


def main():
    """CLI entrypoint dispatching to LINCS or cell-type TWAS mode."""
    argv = sys.argv[1:]
    if argv and argv[0] == "lincs":
        run_lincs_mode(argv[1:])
        return
    if "--lincs" in argv:
        run_lincs_mode([a for a in argv if a != "--lincs"])
        return
    if argv and argv[0] == "cell-type-twas":
        argv = ["--cell-type-twas"] + argv[1:]
    elif argv and (argv[0] in ALL_STEPS or argv[0] in STEP_ALIASES):
        argv = ["--steps", argv[0]] + argv[1:]
    run_twas_mode(argv)


if __name__ == "__main__":
    main()
