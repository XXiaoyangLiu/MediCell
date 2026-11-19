#!/usr/bin/env python3
"""Run an exported bash function across items in parallel with fail-fast semantics."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description="Run exported bash function in parallel")
    p.add_argument("--cores", type=int, default=1)
    p.add_argument("--fn", required=True, help="exported bash function name")
    p.add_argument("items", nargs=argparse.REMAINDER)
    args = p.parse_args()

    items = [x for x in args.items if x != "--"]
    if not items:
        return 0
    cores = max(1, min(args.cores, len(items)))

    def run_one(item):
        proc = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", f"{args.fn} \"$1\"", "_", item],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return item, proc.returncode, proc.stdout

    if cores == 1:
        for item in items:
            _, rc, out = run_one(item)
            if out:
                sys.stdout.write(out)
                sys.stdout.flush()
            if rc != 0:
                return rc
        return 0

    with ThreadPoolExecutor(max_workers=cores) as pool:
        futures = [pool.submit(run_one, item) for item in items]
        for fut in as_completed(futures):
            item, rc, out = fut.result()
            if out:
                sys.stdout.write(out)
                sys.stdout.flush()
            if rc != 0:
                sys.stderr.write(f"Parallel worker failed for item={item} (exit {rc})\n")
                return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
