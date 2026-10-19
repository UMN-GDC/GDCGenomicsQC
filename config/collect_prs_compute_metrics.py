#!/usr/bin/env python3
"""Collect Slurm resource accounting for a manuscript PRS benchmark run."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


FIELDS = [
    "benchmark_stage",
    "job_id",
    "rule",
    "method_group",
    "state",
    "exit_code",
    "elapsed_seconds",
    "elapsed_hours",
    "allocated_cpus",
    "allocated_core_hours",
    "total_cpu",
    "max_rss",
    "max_rss_bytes",
    "requested_memory",
    "start",
    "end",
    "slurm_log",
    "controller_log",
]


def parse_memory(value: str) -> int | None:
    value = value.strip()
    if not value or value in {"Unknown", "0"}:
        return None
    match = re.fullmatch(r"([0-9.]+)([KMGTP]?)", value, flags=re.IGNORECASE)
    if not match:
        return None
    factor = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "P": 1024**5}
    return round(float(match.group(1)) * factor[match.group(2).upper()])


def method_group(rule: str) -> str:
    mapping = {
        "runSingleAncestryPRS": "single_ancestry_bundle",
        "runMultiAncestryPRSCSx": "PRS-CSx",
        "runMultiAncestryCTSLEB": "CT-SLEB",
        "runMultiAncestryPROSPER": "PROSPER",
        "simulateBivariatePhenotypes": "phenotype_simulation",
        "preparePRSInputs": "shared_input_preparation",
        "preparePROSPERTuningTestingSets": "PROSPER_input_preparation",
        "preparePROSPERLassosumParams": "PROSPER_parameter_preparation",
    }
    return mapping.get(rule, rule or "unmapped")


def read_job_ids(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"benchmark_stage", "job_id", "controller_log"}
    if not rows:
        return []
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Job-id file missing columns: {', '.join(sorted(missing))}")
    seen: set[str] = set()
    return [row for row in rows if row["job_id"] and not (row["job_id"] in seen or seen.add(row["job_id"]))]


def sacct(job_id: str) -> list[dict[str, str]]:
    command = [
        "sacct", "-P", "-n", "-j", job_id,
        "--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,ReqMem,Start,End",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        return []
    columns = ["job_id_raw", "state", "exit_code", "elapsed_raw", "alloc_cpus", "total_cpu", "max_rss", "requested_memory", "start", "end"]
    rows = []
    for line in completed.stdout.splitlines():
        fields = line.split("|")
        if len(fields) < len(columns):
            continue
        rows.append(dict(zip(columns, fields)))
    return rows


def rule_log(workflow_dir: Path, job_id: str) -> tuple[str, str]:
    matches = list((workflow_dir / ".snakemake" / "slurm_logs").glob(f"rule_*/*{job_id}.log"))
    if not matches:
        return "", ""
    path = matches[0]
    rule = path.parent.name.removeprefix("rule_")
    return rule, str(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--job-id-file", required=True, type=Path)
    args = parser.parse_args()

    job_rows = read_job_ids(args.job_id_file)
    output_rows: list[dict[str, str | int | float]] = []
    for job in job_rows:
        accounting = sacct(job["job_id"])
        main_row = next((row for row in accounting if row["job_id_raw"] == job["job_id"]), None)
        if main_row is None:
            main_row = {key: "" for key in (accounting[0] if accounting else {"job_id_raw": ""})}
        max_rss = max((parse_memory(row.get("max_rss", "")) or 0 for row in accounting), default=0)
        rule, log_path = rule_log(args.workflow_dir, job["job_id"])
        elapsed = int(main_row.get("elapsed_raw", "0") or 0) if str(main_row.get("elapsed_raw", "")).isdigit() else 0
        cpus = int(main_row.get("alloc_cpus", "0") or 0) if str(main_row.get("alloc_cpus", "")).isdigit() else 0
        output_rows.append({
            "benchmark_stage": job["benchmark_stage"],
            "job_id": job["job_id"],
            "rule": rule,
            "method_group": method_group(rule),
            "state": main_row.get("state", "accounting_unavailable"),
            "exit_code": main_row.get("exit_code", ""),
            "elapsed_seconds": elapsed or "",
            "elapsed_hours": round(elapsed / 3600, 4) if elapsed else "",
            "allocated_cpus": cpus or "",
            "allocated_core_hours": round(elapsed * cpus / 3600, 4) if elapsed and cpus else "",
            "total_cpu": main_row.get("total_cpu", ""),
            "max_rss": max((row.get("max_rss", "") for row in accounting), key=lambda value: parse_memory(value) or 0, default=""),
            "max_rss_bytes": max_rss or "",
            "requested_memory": main_row.get("requested_memory", ""),
            "start": main_row.get("start", ""),
            "end": main_row.get("end", ""),
            "slurm_log": log_path,
            "controller_log": job["controller_log"],
        })

    destination = args.output_dir / "benchmark_compute_metrics.tsv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {destination} with {len(output_rows)} Slurm jobs.")


if __name__ == "__main__":
    main()
