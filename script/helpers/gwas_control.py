#!/usr/bin/env python3
"""Run METAL once for either validated standard studies or the isolated SLE adapter.
gwas_run.json is owned by the driver; diagnostics live in metal.log and stderr."""
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys

LEGACY = ("sle-legacy-pbmc", "sle-legacy-kidney")
DIAGNOSTICS = ["183 BETA is column 26 zscore_abs, not column 25 beta", "183 PBMC SE is column 28 signed z; kidney SE is column 8", "096 default awk whitespace consumes empty fields", "ALLELE A2 A1 reverses effect allele declaration", "SAMPLESIZE N is an invalid METAL command", "OR/P uses historical 1-P/2 formula", "MC_CHROMS is intentionally ignored; chromosomes 1-22", "METAL rounded 4-digit values feed legacy filtering"]

def is_known_legacy_warning(line):
    return bool(re.fullmatch(
        r"## WARNING: An additional [0-9]+ rows with duplicate marker names were ignored",
        line.strip())) or any(t in line.lower() for t in (
            "command you issued could not be processed", "samplesize",
            "invalid standard error", "invalid value for standard error",
            "bad alleles", "bad allele", "invalid effect", "invalid pvalue",
            "invalid frequency", "missing value for marker", "invalid marker",
            "allele mismatch", "duplicated in input"))


def validate_result(target, legacy=False):
    required = {"MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value",
                "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal"}
    count = 0
    seen = set()
    with target.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Unexpected METAL output schema")
        for row in reader:
            marker = row["MarkerName"]
            if not marker or marker in seen:
                raise ValueError("Missing/duplicate METAL marker")
            seen.add(marker)
            values = {key: float(row[key]) for key in (
                "Effect", "StdErr", "P-value", "HetISq", "HetChiSq", "HetDf", "HetPVal")}
            if (not all(math.isfinite(v) for v in values.values())
                    or values["StdErr"] <= 0 or not 0 <= values["P-value"] <= 1
                    or not 0 <= values["HetISq"] <= 100 or values["HetChiSq"] < 0
                    or values["HetDf"] < 0 or not values["HetDf"].is_integer()
                    or not 0 <= values["HetPVal"] <= 1
                    or not row["Allele1"] or not row["Allele2"]
                    or (not legacy and (row["Allele1"].upper() not in ("A", "C", "G", "T")
                        or row["Allele2"].upper() not in ("A", "C", "G", "T")
                        or row["Allele1"].upper() == row["Allele2"].upper()))):
                raise ValueError("Invalid numeric/allele METAL result: " + marker)
            count += 1
    if not count:
        raise ValueError("No METAL result variants")
    return count


def validate_normalized(path):
    seen = set()
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != ["SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"]:
            raise ValueError("Invalid normalized study schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Invalid normalized study row width")
            key = row["SNP"]
            if key in seen or not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2])_[1-9][0-9]*_[ACGT]_[ACGT]", key):
                raise ValueError("Invalid/duplicate normalized variant key")
            seen.add(key)
            parts = key.split("_")
            if row["A1"] != parts[3] or row["A2"] != parts[2] or parts[2] == parts[3]:
                raise ValueError("Invalid normalized effect allele/key mapping")
            if not math.isfinite(float(row["BETA"])) or not math.isfinite(float(row["SE"])) or float(row["SE"]) <= 0:
                raise ValueError("Invalid normalized BETA/SE")
            for column in ("P", "FREQ", "N"):
                if row[column] == "NA":
                    continue
                value = float(row[column])
                if (not math.isfinite(value) or (column == "N" and value <= 0)
                        or (column != "N" and not 0 <= value <= 1)):
                    raise ValueError("Invalid normalized " + column)
    if not seen:
        raise ValueError("Empty normalized study")


def config():
    e = os.environ
    profile, manifest = e.get("MC_GWAS_PROFILE", ""), e.get("MC_GWAS_MANIFEST", "")
    if bool(profile) == bool(manifest) or (profile and profile not in LEGACY):
        raise ValueError("Choose exactly one supported MC_GWAS_PROFILE or MC_GWAS_MANIFEST")
    out = Path(e["MC_GWAS_OUT"]).resolve()
    prefix = e["MC_GWAS_PREFIX"]
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", prefix):
        raise ValueError("Unsafe MC_GWAS_PREFIX")
    if profile and prefix != "SLE_metal_096_917_183":
        raise ValueError("Legacy prefix must be SLE_metal_096_917_183")
    policy = e.get("MC_GWAS_FILTER_POLICY") or ("legacy" if profile else "none")
    if policy not in ("none", "legacy") or (profile and policy != "legacy"):
        raise ValueError("Invalid GWAS filter policy")
    if Path(e["MC_METAL_TBL"]).resolve() != out / (prefix + "1.tbl"):
        raise ValueError("MC_METAL_TBL must equal <out>/<prefix>1.tbl")
    return out, prefix, profile

def meta():
    out, prefix, profile = config()
    if profile:
        print("WARNING: LEGACY REPRODUCTION ONLY; known statistical errors intentionally retained: " + "; ".join(DIAGNOSTICS), flush=True)
        subprocess.run(["bash", str(Path(__file__).with_name("legacy_sle_metal.sh"))], check=True)
    else:
        with (out / "gwas_studies.tsv").open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames != ["study_id", "path"]:
                raise ValueError("Invalid normalized study manifest header")
            studies = list(reader)
        if (not studies or len({s["study_id"] for s in studies}) != len(studies)
                or any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", s["study_id"])
                       for s in studies)):
            raise ValueError("No studies or invalid/duplicate normalized study IDs")
        # P/FREQ/N are not used by SCHEME STDERR; declaring them only invites METAL warnings.
        commands = ["SCHEME STDERR", "SEPARATOR TAB", "MARKER SNP", "ALLELE A1 A2", "EFFECT BETA", "STDERR SE"]
        for s in studies:
            p = Path(s["path"])
            # METAL command parser has no safe filename quoting contract; reject whitespace, not shell-eval it.
            if (not p.is_absolute() or not p.is_file() or re.search(r"\s", str(p))
                    or p.resolve().parent != out
                    or p.name != s["study_id"] + "_normalized.tsv"):
                raise ValueError("Normalized study paths must exist and contain no whitespace")
            validate_normalized(p)
            commands.append("PROCESS " + str(p))
        commands += ["OUTFILE " + prefix + " .tbl", "ANALYZE HETEROGENEITY", "QUIT"]
        (out / "metal_script.txt").write_text("\n".join(commands) + "\n")
    target = out / (prefix + "1.tbl")
    # Fail rather than accept stale output or let METAL silently switch to ...2.tbl.
    if target.exists() or target.with_suffix(".tbl.info").exists():
        raise ValueError("METAL target already exists; use a fresh output directory or explicitly remove own previous products")
    binary = os.environ.get("METALBIN") or os.environ.get("MC_METAL_BIN")
    if not binary:
        raise ValueError("METALBIN or MC_METAL_BIN required")
    log = out / "metal.log"
    with log.open("w") as f:
        result = subprocess.run([binary, "metal_script.txt"], cwd=out, stdout=f, stderr=subprocess.STDOUT)
    text = log.read_text()
    problems = []
    for line in text.splitlines():
        if re.search(r"\b(?:error|warning)\s*[:\-]|failed to|unknown command|unrecognized command", line.lower()):
            # Full diagnostic lines stay in metal.log. Only the known SLE defects are tolerated here.
            expected = profile and is_known_legacy_warning(line)
            if not expected:
                problems.append(line)
    if result.returncode or problems or not target.is_file() or target.stat().st_size == 0:
        raise RuntimeError("METAL failed internal/output validation; inspect " + str(log) + ": " + repr(problems[:8]))
    count = validate_result(target, legacy=bool(profile))
    print("METAL validated: %d variants; %s" % (count, target))

if __name__ == "__main__":
    try:
        meta()
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        sys.exit("GWAS: " + str(exc))
