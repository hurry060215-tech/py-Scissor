"""Explicit developer command to retrieve the complete official LUAD inputs.

Never imported by the library. No mirrors, credentials, retries, or sampling.
The original endpoints are documented in the fixed Scissor tutorial.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import urllib.request

BASE = "https://xialab.s3-us-west-2.amazonaws.com/Duanchen/Scissor_data/"
FILES = (
    "TCGA_LUAD_survival.RData", "TCGA_LUAD_TP53_mutation.RData",
    "scRNA-seq.RData", "TCGA_LUAD_exp1.RData", "TCGA_LUAD_exp2.RData",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    manifest = {
        "source": "https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html",
        "requested_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "sampling": "none", "files": [],
    }
    record = args.destination / "download-manifest.json"
    for name in FILES:
        path = args.destination / name
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}; verify existing data separately")
        part = path.with_suffix(path.suffix + ".part")
        item = {"name": name, "url": BASE + name}
        start = time.monotonic()
        try:
            with urllib.request.urlopen(BASE + name, timeout=60) as response:
                item.update(status=response.status,
                            expected_bytes=int(response.headers["Content-Length"]),
                            etag=response.headers.get("ETag"),
                            last_modified=response.headers.get("Last-Modified"))
                if item["expected_bytes"] > 500_000_000:
                    raise RuntimeError("Unexpectedly large official file; inspect before retrieval")
                digest = hashlib.sha256()
                size = 0
                with part.open("xb") as out:
                    while chunk := response.read(1024 * 1024):
                        if time.monotonic() - start > 600:
                            raise TimeoutError("Per-file 600 second total time budget exceeded")
                        size += len(chunk)
                        if size > item["expected_bytes"]:
                            raise RuntimeError("Response exceeds announced Content-Length")
                        out.write(chunk)
                        digest.update(chunk)
                if size != item["expected_bytes"]:
                    raise RuntimeError("Incomplete response")
                part.replace(path)
                item.update(bytes=size, sha256=digest.hexdigest(), complete=True)
        except Exception as error:
            item.update(complete=False, error_type=type(error).__name__, error=str(error))
            raise
        finally:
            item["elapsed_seconds"] = round(time.monotonic() - start, 3)
            manifest["files"].append(item)
            record.write_text(json.dumps(manifest, indent=2) + "\n")
            print(json.dumps(item), flush=True)


if __name__ == "__main__":
    main()
