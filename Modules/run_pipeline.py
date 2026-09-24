"""Orchestrator script that triggers Snakemake."""

import argparse
import os
import subprocess
import sys
from logging import getLogger, basicConfig, INFO

basicConfig(level=INFO, format="[%(levelname)s] %(message)s")
logger = getLogger("Orchestrator")


def main():
    parser = argparse.ArgumentParser(
        description="Run the ReCoDe Pipeline via Snakemake"
    )
    parser.add_argument("config_file", help="Path to the TOML config file.")
    parser.add_argument("--cores", type=int, default=32, help="Number of cores to use.")
    parser.add_argument(
        "--dry-run", action="store_true", help="Print the DAG without running anything."
    )
    args = parser.parse_args()

    os.environ["RECODE_CONFIG"] = args.config_file

    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        "Modules/Snakefile",
        "--cores",
        str(args.cores),
        "--keep-going",
        "--printshellcmds",
    ]

    if args.dry_run:
        cmd.append("--dry-run")

    logger.info(f"Starting pipeline orchestration with {args.cores} cores...")
    logger.info(f"Command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        if args.dry_run:
            logger.info("Dry run successful! The DAG is perfectly connected.")
        else:
            logger.info("Pipeline completed successfully!")
    except subprocess.CalledProcessError:
        logger.error(
            "Pipeline failed! Check the logs above for the specific module error."
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
