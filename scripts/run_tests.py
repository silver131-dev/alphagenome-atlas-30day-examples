#!/usr/bin/env python3
"""Run self-contained Day02-Day30 unittest modules with their local fixture paths.
No pip installation, SDK import, API credentials, or network needed.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DIRS=sorted([p for p in (ROOT/'examples').iterdir() if p.is_dir()])

def run_one(p:Path, verbose:bool=False)->int:
    pattern='test_*.py'
    cmd=[sys.executable,'-B','-m','unittest','discover','-s','.','-p',pattern]
    if verbose:cmd.append('-v')
    result=subprocess.run(cmd,cwd=p,text=True,capture_output=True,timeout=240)
    print(f"{p.relative_to(ROOT)}: {'PASS' if result.returncode==0 else 'FAIL'}")
    if verbose or result.returncode:
        print(result.stdout[-3000:])
        print(result.stderr[-4500:])
    else:
        for line in result.stderr.splitlines():
            if line.startswith(('Ran ','OK')):print('  '+line)
    return result.returncode

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--only',help='day03, variant_pipeline, day29 or historical_readers')
    ap.add_argument('--verbose',action='store_true')
    args=ap.parse_args()
    folders=[p for p in DIRS if not args.only or p.name==args.only]
    if not folders:ap.error('No matching folder')
    fail=sum(run_one(p,args.verbose)!=0 for p in folders)
    print(f"DIRECTORIES={len(folders)} FAILURES={fail} REAL_API_CALLS_THIS_RUN=0")
    return 1 if fail else 0
if __name__=='__main__':raise SystemExit(main())
