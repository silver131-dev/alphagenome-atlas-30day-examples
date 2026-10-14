#!/usr/bin/env python3
"""Day02 OPTIONAL synthetic biosample data-contract lesson; never calls APIs."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def validate_record(data: dict) -> dict:
    if not isinstance(data, dict) or data.get("data_kind") != "SYNTHETIC_EDUCATIONAL_DEMO":
        raise ValueError("Only explicitly synthetic teaching data allowed")
    if data.get("genome_build") != "TOY_NOT_GRCH38" or not str(data.get("variant", "")).startswith("chrToy:"):
        raise ValueError("Real genome identifiers are forbidden")
    if data.get("model_prediction_performed") is not False or data.get("real_api_calls_this_run") != 0:
        raise ValueError("No model/API data accepted")
    if data.get("assay") != "splicing_proxy_delta":
        raise ValueError("Assay must match the teaching data contract")
    rows=data.get("rows")
    if not isinstance(rows,list) or len(rows)!=3:
        raise ValueError("Three fictional biosamples required")
    names=set()
    result=[]
    for row in rows:
        if not isinstance(row,dict):
            raise ValueError("Invalid row")
        name=row.get("biosample")
        if name not in {"Neuron-like","Blood-like","Liver-like"} or name in names:
            raise ValueError("Unexpected or repeated biosample")
        names.add(name)
        if row.get("experimental_observation")!="NOT_PERFORMED":
            raise ValueError("No observed experiments in the demo")
        score=row.get("score")
        status=row.get("score_status")
        if score is None:
            if status!="MISSING_SCORE":
                raise ValueError("Null must remain missing")
            label="MISSING_NOT_ZERO"
        else:
            if isinstance(score,bool) or not isinstance(score,(int,float)) or not math.isfinite(score) or not -1<=score<=1 or status!="SYNTHETIC_NUMBER":
                raise ValueError("Invalid teaching score")
            label="SYNTHETIC_ZERO" if score==0 else "SYNTHETIC_NONZERO"
        result.append({"biosample":name,"score":score,"score_status":status,"interpretation":label,"experimental_observation":"NOT_PERFORMED"})
    if names!={"Neuron-like","Blood-like","Liver-like"}:
        raise ValueError("Missing biosample")
    return {"status":"PASS","mode":"OFFLINE_ONLY","data_kind":"SYNTHETIC_EDUCATIONAL_DEMO",
      "genome_build":"TOY_NOT_GRCH38","variant":data["variant"],"assay":data["assay"],
      "rows":result,"null_count":sum(x["score"] is None for x in result),
      "model_prediction_performed":False,"experimental_validation_performed":False,
      "clinical_causality":"NOT_ESTABLISHED","manual_review":"REQUIRED","real_api_calls_this_run":0}


def run(fixture: Path, output_dir: Path) -> dict:
    # Read only user-supplied LOCAL synthetic JSON; no network functionality exists.
    raw=fixture.read_bytes()
    data=json.loads(raw)
    packet=validate_record(data)
    output_dir.mkdir(parents=True,exist_ok=True)
    result_path=output_dir/'day02_biosample_card.json'
    csv_path=output_dir/'day02_biosample_card.csv'
    result_path.write_text(json.dumps(packet,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    with csv_path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(packet['rows'][0]))
        w.writeheader();w.writerows(packet['rows'])
    return {"status":"PASS","mode":"OFFLINE_ONLY","biosample_count":len(packet['rows']),
      "null_count":packet['null_count'],"source_sha256":hashlib.sha256(raw).hexdigest(),
      "json_file":str(result_path),"csv_file":str(csv_path),"real_api_calls_this_run":0}


def main():
    p=argparse.ArgumentParser(description="Day 02 optional synthetic biosample card")
    p.add_argument('--mode',required=True,choices=['offline'])
    p.add_argument('--fixture',default='day02_synthetic_biosample_card.json')
    p.add_argument('--output-dir',default='day02_output')
    x=p.parse_args()
    print(json.dumps(run(Path(x.fixture),Path(x.output_dir)),ensure_ascii=False,sort_keys=True))

if __name__=='__main__':
    main()
