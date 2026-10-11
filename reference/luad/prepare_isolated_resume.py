"""GPL-3.0-or-later. Snapshot completed Cox nodes into a new, isolated root.

Explicitly preserves the fact that the old writer's lifecycle is unknown.
Copies only full terminal nodes, never partially completed model paths.
Independent R certificates must pass before computing any missing node.
"""
import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("destination", type=Path)
parser.add_argument("--original-session-id", type=int, default=None,
                    help="Optional local execution-channel identifier for private provenance")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
base = root / "validation-data/luad/r-oracle"
source = base / "corrected-tight"
destination = args.destination.resolve()
assert not destination.exists(), "Isolation requires a previously nonexistent output root"
assert destination.is_relative_to(base) and destination != source
manifest_path = base / "tight-cox-inputs/manifest.json"
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest_sha = sha(manifest_path)
assert manifest_sha == "a8cd4b89308595d83fdcb02e9edbd5408d4d430f9c8de17a12ee607e76bb465d"
manifest = json.loads(manifest_path.read_text())
source_log = base / "corrected-tight-cox.log"
log_hash = sha(source_log)
completed = {"0.05": ["full", "fold-1", "fold-2", "fold-3", "fold-4", "fold-5", "refit"],
             "0.2": ["full", "fold-1", "fold-2"]}
nodes, files = [], {}
for alpha, roles in completed.items():
    directory = source / "cox" / f"alpha-{alpha}"
    result_path = directory / "result.json"
    if alpha == "0.2":
        assert not result_path.exists(), "Old alpha .2 advanced beyond this fixed recovery boundary"
    result = json.loads(result_path.read_text()) if result_path.exists() else None
    if result is not None:
        assert result["metadata"]["family"] == "cox" and result["metadata"]["thresh"] == 1e-10
        assert all(result["valid"])
    for role in roles:
        checkpoint = directory / f"{role}-checkpoint.json"
        if checkpoint.exists():
            saved = json.loads(checkpoint.read_text())
            json_path = checkpoint
        else:
            assert result is not None
            saved = result[role] if role in ("full", "refit") else result["folds"][int(role[5:]) - 1]
            json_path = result_path
        assert saved["coefficient_file"] == f"{role}-coefficients.f64"
        shape = saved["coefficient_shape"]
        assert shape == [4102, 1 if role == "refit" else 20]
        flag = saved["converged"]
        flags = flag if isinstance(flag, list) else [flag]
        assert len(flags) == shape[1] and all(flags)
        residual = saved["kkt_residual"]
        residuals = residual if isinstance(residual, list) else [residual]
        assert max(residuals) <= 1e-10
        binary = directory / saved["coefficient_file"]
        assert binary.stat().st_size == np.prod(shape) * 8
        assert np.isfinite(np.fromfile(binary, dtype="<f8")).all()
        node_files = [json_path, binary]
        if result is not None:
            node_files.append(result_path)
        for path in node_files:
            name = str(path.relative_to(source))
            files[name] = {"source": str(path), "sha256": sha(path), "bytes": path.stat().st_size}
        nodes.append({"alpha": float(alpha), "role": role, "coefficient_shape": shape,
                      "independent_R_recertification_required": True})
assert len(nodes) == 10 and sha(source_log) == log_hash
destination.mkdir()
for name, record in files.items():
    src = Path(record["source"])
    dst = destination / name
    dst.parent.mkdir(parents=True, exist_ok=True)
    assert sha(src) == record["sha256"]
    shutil.copy2(src, dst)  # Deliberately real copies, never links to an unknown writer.
    assert sha(dst) == sha(src) == record["sha256"]
assert sha(source_log) == log_hash, "Old log advanced during snapshot: do not schedule new fits"
snapshot = {"created_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "copied_hash_verified_pending_independent_R_certificates",
            "source_root": str(source), "output_root": str(destination),
            "original_writer_lifecycle": "unknown_after_executor_key_change",
            "original_session_id": args.original_session_id, "input_manifest_sha256": manifest_sha,
            "source_log": str(source_log), "source_log_sha256": log_hash,
            "source_log_mtime_ns": source_log.stat().st_mtime_ns,
            "completed_nodes": nodes, "files": files,
            "remaining_nodes": {"0.2": ["fold-3", "fold-4", "fold-5", "refit"],
                                "0.6": ["full", "fold-1", "fold-2", "fold-3", "fold-4", "fold-5", "refit"]},
            "old_output_modified": False, "copies_are_independent_regular_files": True}
(destination / "isolation-snapshot.json").write_text(json.dumps(snapshot, indent=2) + "\n")
print(json.dumps({"output": str(destination), "copied_nodes": len(nodes), "copied_files": len(files),
                  "snapshot_sha256": sha(destination / "isolation-snapshot.json")}))
