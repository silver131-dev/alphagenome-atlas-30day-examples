#!/usr/bin/env python3
"""Conservative, offline public-candidate preflight. Not a security certification."""
from __future__ import annotations
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FORBIDDEN_FILES={'.env','.env.local','.env.production','credentials.json','service-account.json','id_rsa','id_ed25519'}
PATTERNS={
    'google_api_key':r'AIza[0-9A-Za-z_-]{20,}',
    'openai_key':r'\bsk-[A-Za-z0-9_-]{20,}',
    'github_token':r'\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{20,}',
    'aws_access_key':r'\bAKIA[A-Z0-9]{16}\b',
    'private_key_pem':r'-----BEGIN(?: RSA| EC| OPENSSH)? PRIVATE KEY-----',
    'jwt':r'\beyJ[A-Za-z0-9_-]{25,}\.eyJ[A-Za-z0-9_-]{20,}\.\w+',
}

def main()->int:
    findings=[];checked=0
    for f in ROOT.rglob('*'):
        if not f.is_file():continue
        if '.git' in f.parts or '__pycache__' in f.parts:continue
        checked+=1
        if f.name in FORBIDDEN_FILES or f.suffix.lower() in {'.pem','.p12','.pfx','.key'}:
            findings.append(f'{f.relative_to(ROOT)}: restricted filename')
        try:s=f.read_text(encoding='utf-8')
        except UnicodeDecodeError:continue
        for name,pattern in PATTERNS.items():
            if re.search(pattern,s):findings.append(f'{f.relative_to(ROOT)}: possible {name}')
    print(f'SCANNED_FILES={checked} SUSPECT_FINDINGS={len(findings)}')
    for line in findings:print(line)
    print('LIMITED_REGEX_SCAN_ONLY=TRUE (not proof of no PHI/PII/secrets)')
    return 1 if findings else 0
if __name__=='__main__':raise SystemExit(main())
