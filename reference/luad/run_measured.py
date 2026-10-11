"""Run a reference command and record actual exit/elapsed/child peak RSS."""
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

destination, *command = sys.argv[1:]
start = time.monotonic()
result = subprocess.run(command, check=False)
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
summary = {"command": command, "exit_code": result.returncode,
           "elapsed_seconds": time.monotonic() - start,
           "child_peak_rss_bytes": usage.ru_maxrss * 1024,
           "user_cpu_seconds": usage.ru_utime, "system_cpu_seconds": usage.ru_stime}
Path(destination).write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
sys.exit(result.returncode)
