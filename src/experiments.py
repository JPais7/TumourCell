"""Reproducible experiment manifests."""
from __future__ import annotations
import hashlib,json,platform,subprocess
from datetime import datetime,timezone
from pathlib import Path
def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()
def write_manifest(experiment_id,*,config,dataset_manifest,inputs,outputs,model_version,random_seed,root=Path("results/experiments")):
    out=root/experiment_id; out.mkdir(parents=True,exist_ok=False)
    try: commit=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()
    except Exception as e: raise RuntimeError("git commit unavailable") from e
    payload={"experiment_id":experiment_id,"config":config,"git_commit":commit,"dataset_manifest":str(dataset_manifest),
      "input_hashes":{str(p):sha256(p) for p in inputs},"model_version":model_version,"random_seed":random_seed,
      "software_environment":{"python":platform.python_version(),"platform":platform.platform()},
      "output_hashes":{str(p):sha256(p) for p in outputs},"timestamp":datetime.now(timezone.utc).isoformat()}
    (out/"manifest.json").write_text(json.dumps(payload,indent=2)+"\n"); return payload
