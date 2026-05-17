"""Runtime tracking module."""

import time
import pandas as pd
from pathlib import Path
from contextlib import contextmanager
from logging import getLogger

logger = getLogger(__name__)


class RuntimeTracker:
    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.records = []

        # Ensure the directory exists
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

    @contextmanager
    def measure(self, process_name: str):
        """Context manager to measure execution time of a block of code."""
        logger.info(f"[START] Tracking runtime for: {process_name}")
        start_time = time.time()
        status = "Completed"

        try:
            yield
        except Exception:
            status = "Failed"
            raise
        finally:
            end_time = time.time()
            duration_sec = end_time - start_time
            duration_min = duration_sec / 60.0

            # Add record to list
            self.records.append(
                {
                    "Process": process_name,
                    "Status": status,
                    "Duration (Seconds)": round(duration_sec, 2),
                    "Duration (Minutes)": round(duration_min, 2),
                }
            )

            # Save incrementally
            self.save()

            # Log the result
            if status == "Completed":
                logger.info(
                    f"[DONE] {process_name} finished in {duration_sec:.2f}s ({duration_min:.2f}m)"
                )
            else:
                logger.error(
                    f"[FAILED] {process_name} crashed after {duration_sec:.2f}s ({duration_min:.2f}m)"
                )

    def save(self):
        """Saves the runtime records to a CSV file."""
        df = pd.DataFrame(self.records)
        df.to_csv(self.output_path, index=False)
