"""Day 07: small, hand-crafted CNN + scalar attention educational analogy.

No machine learning is trained. No AlphaGenome weights, SDK, remote network,
patient information, model output, or official internal architecture are used.
Python 3.12+ standard library only. Input must be a controlled synthetic fixture.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
import math
import re
from pathlib import Path
from typing import Any

EVIDENCE = "SYNTHETIC_EDUCATIONAL_DEMO"
OUTPUT_FILES = (
    "day07_multiscale_report.json",
    "day07_token_attention.csv",
    "day07_multiscale_map.svg",
    "day07_multiscale_dashboard.html",
)
MANIFEST_NAME = "day07_sha256_manifest.json"
DNA = set("ACGT")
EXPECTED_TOP_LEVEL = {
    "evidence_class", "assembly", "chrom", "sequence_start_1based",
    "sequence", "bin_size", "local_motif", "distal_motif",
    "local_motif_start_1based", "distal_motif_start_1based", "variants",
}
EXPECTED_VARIANT_FIELDS = {
    "label", "variant_id", "position_1based", "ref", "alt", "source",
}
LABELS = ("DISTAL_CHANGE", "LOCAL_CHANGE")


def _integer(v: Any) -> bool:
    return type(v) is int


def load_fixture(path: str | Path) -> tuple[dict[str, Any], str]:
    """Load strictly controlled 64-bp demonstration, never a patient sequence."""
    target = Path(path)
    raw = target.read_bytes()
    data = json.loads(raw)
    if type(data) is not dict or set(data) != EXPECTED_TOP_LEVEL:
        raise ValueError("Expected the Day07 synthetic fixture schema")
    if data["evidence_class"] != EVIDENCE or data["assembly"] != "TOY_NOT_GRCH38":
        raise ValueError("Only explicitly toy data are permitted")
    if data["chrom"] != "chrToy" or not _integer(data["sequence_start_1based"]):
        raise ValueError("Expected chrToy and an integer sequence start")
    if data["sequence_start_1based"] <= 0:
        raise ValueError("Sequence start must be positive")
    sequence = data["sequence"]
    if type(sequence) is not str or len(sequence) != 64 or not set(sequence) <= DNA:
        raise ValueError("Expected precisely 64 upper-case synthetic DNA bases")
    if not _integer(data["bin_size"]) or data["bin_size"] != 8:
        raise ValueError("The toy design uses eight 8-bp bins")
    for k in ("local_motif", "distal_motif"):
        m = data[k]
        if type(m) is not str or len(m) != 3 or not set(m) <= DNA:
            raise ValueError("Motif must contain three upper-case DNA letters")
    for key, motif in (("local_motif_start_1based", "local_motif"),
                       ("distal_motif_start_1based", "distal_motif")):
        pos = data[key]
        i = pos - data["sequence_start_1based"] if _integer(pos) else -1
        if i < 0 or i+3 > 64 or sequence[i:i+3] != data[motif]:
            raise ValueError("Annotated motif must exactly match toy DNA")
        if find_exact_hits(sequence, data[motif]) != [i]:
            raise ValueError("Only a single full hit of each motif is supported")
    if data["distal_motif_start_1based"]-data["local_motif_start_1based"] != 40:
        raise ValueError("Toy long-range distance must remain 40 bp")
    var = data["variants"]
    if type(var) is not list or len(var) != 2:
        raise ValueError("Two predefined variants are required")
    if [v.get("label") if type(v) is dict else None for v in var] != list(LABELS):
        raise ValueError("Expected distal and local alterations in that order")
    for v in var:
        if type(v) is not dict or set(v) != EXPECTED_VARIANT_FIELDS:
            raise ValueError("Unexpected variant fields")
        if v["source"] != "SYNTHETIC_FIXTURE":
            raise ValueError("Input is not a synthetic variant")
        pos = v["position_1based"]
        if not _integer(pos) or not data["sequence_start_1based"] <= pos < data["sequence_start_1based"]+64:
            raise ValueError("Variant lies outside toy sequence")
        i = pos-data["sequence_start_1based"]
        if type(v["ref"]) is not str or type(v["alt"]) is not str or len(v["ref"]) != 1 or len(v["alt"]) != 1:
            raise ValueError("SNV must be single-base")
        if v["ref"] not in DNA or v["alt"] not in DNA or v["ref"] == v["alt"]:
            raise ValueError("Expected a proper synthetic SNV")
        if sequence[i] != v["ref"] or v["variant_id"] != f"chrToy:{pos}:{v['ref']}>{v['alt']}":
            raise ValueError("REF or identifier does not match toy sequence")
        expected_index = (data["distal_motif_start_1based"] if v["label"] == "DISTAL_CHANGE"
                          else data["local_motif_start_1based"]+1)
        if pos != expected_index:
            raise ValueError("SNV must alter the designated motif")
        mutated = apply_snv(sequence, pos, v["ref"], v["alt"], data["sequence_start_1based"])
        lost = data["distal_motif"] if v["label"] == "DISTAL_CHANGE" else data["local_motif"]
        if find_exact_hits(mutated, lost):
            raise ValueError("Toy mutation must erase its target motif")
        other = (data["local_motif"] if v["label"] == "DISTAL_CHANGE" else data["distal_motif"])
        if find_exact_hits(mutated, other) != find_exact_hits(sequence, other):
            raise ValueError("Target mutation must not create or erase the other teaching motif")
    return data, hashlib.sha256(raw).hexdigest()


def find_exact_hits(sequence: str, motif: str) -> list[int]:
    """Full-match positions in zero-based local sequence coordinates."""
    if not motif or not set(motif) <= DNA or not set(sequence) <= DNA:
        raise ValueError("Expected nonempty DNA motif and DNA sequence")
    return [i for i in range(len(sequence)-len(motif)+1)
            if sequence[i:i+len(motif)] == motif]


def fixed_kernel_scan(sequence: str, motif: str) -> list[float]:
    """Untrained, one-hot-like width-3 sliding comparison (not a learned CNN)."""
    if type(sequence) is not str or type(motif) is not str or not motif or not set(motif) <= DNA or not set(sequence) <= DNA:
        raise ValueError("Only A/C/G/T upper-case sequences are allowed")
    return [sum(sequence[i+k] == motif[k] for k in range(len(motif)))/len(motif)
            for i in range(len(sequence)-len(motif)+1)]


def pool_full_hits(sequence: str, motif: str, bin_size: int) -> tuple[list[int], list[int]]:
    if not _integer(bin_size) or bin_size <= 0 or len(sequence) % bin_size:
        raise ValueError("Invalid bin size")
    kernel_scores = fixed_kernel_scan(sequence, motif)
    hits = [i for i, score in enumerate(kernel_scores) if score == 1.0]
    profile = [0]*(len(sequence)//bin_size)
    for i in hits:
        profile[i//bin_size] = 1
    return profile, hits


def softmax(values: list[float]) -> list[float]:
    if not values or any(not math.isfinite(float(x)) for x in values):
        raise ValueError("Attention logits must be finite and nonempty")
    center = max(values)
    exps = [math.exp(x-center) for x in values]
    total = sum(exps)
    return [x/total for x in exps]


def scalar_attention(local_tokens: list[int], distal_tokens: list[int], query_bin: int) -> tuple[list[float], float]:
    """Hand-set scalar Q=3*local, K=distal, V=distal; no learned weights."""
    if not local_tokens or len(local_tokens) != len(distal_tokens) or not 0 <= query_bin < len(local_tokens):
        raise ValueError("Invalid attention token inputs")
    if any(v not in (0,1) or type(v) is not int for v in local_tokens+distal_tokens):
        raise ValueError("Expected binary toy pooled motif hits")
    q = 3*local_tokens[query_bin]
    logits = [q*k for k in distal_tokens]
    weights = softmax(logits)
    context = sum(w*v for w, v in zip(weights, distal_tokens))
    return weights, context


def apply_snv(sequence: str, pos: int, ref: str, alt: str, start: int) -> str:
    if not _integer(pos) or not _integer(start) or not (0 <= pos-start < len(sequence)):
        raise ValueError("Variant outside toy sequence")
    if sequence[pos-start] != ref or ref == alt or alt not in DNA:
        raise ValueError("Invalid SNV REF/ALT")
    return sequence[:pos-start] + alt + sequence[pos-start+1:]


def analyze_sequence(sequence: str, fixture: dict[str, Any]) -> dict[str, Any]:
    local, local_hits = pool_full_hits(sequence, fixture["local_motif"], fixture["bin_size"])
    distal, distal_hits = pool_full_hits(sequence, fixture["distal_motif"], fixture["bin_size"])
    query_bin = (fixture["local_motif_start_1based"]-fixture["sequence_start_1based"])//fixture["bin_size"]
    weights, context = scalar_attention(local, distal, query_bin)
    # Toy U-Net analogy: a precise local hit (skip) gated by a coarse attention context.
    # Neither branch is a trained biological predictor.
    signal = local[query_bin]*(0.5+0.5*context)
    decoder_per_base = [0.0]*len(sequence)
    for index in local_hits:
        bin_id = index//fixture["bin_size"]
        _, ctx = scalar_attention(local, distal, bin_id)
        decoder_per_base[index] = round(0.5+0.5*ctx, 6)
    return {
        "local_kernel_full_hits_1based": [fixture["sequence_start_1based"]+i for i in local_hits],
        "distal_kernel_full_hits_1based": [fixture["sequence_start_1based"]+i for i in distal_hits],
        "local_pooled_8bp": local,
        "distal_pooled_8bp": distal,
        "query_token_bin": query_bin,
        "attention_from_local_bin": [round(w, 8) for w in weights],
        "long_range_context_proxy": round(context, 8),
        "toy_fused_output": round(signal, 6),
        "decoder_proxy_1bp_nonzero": [
            {"position_1based": fixture["sequence_start_1based"]+i, "value": v}
            for i, v in enumerate(decoder_per_base) if v > 0
        ],
    }


def build_report(fixture: dict[str, Any], sha256: str, fixture_name: str) -> dict[str, Any]:
    base = analyze_sequence(fixture["sequence"], fixture)
    cases: list[dict[str, Any]] = [{"label":"REFERENCE","variant_id":None,"local_context":base}]
    for variant in fixture["variants"]:
        alt_sequence=apply_snv(fixture["sequence"],variant["position_1based"],variant["ref"],variant["alt"],fixture["sequence_start_1based"])
        case = analyze_sequence(alt_sequence,fixture)
        cases.append({"label":variant["label"],"variant_id":variant["variant_id"],"local_context":case,
                      "delta_alt_minus_ref_toy_output":round(case["toy_fused_output"]-base["toy_fused_output"],6),
                      "local_cnn_unchanged":case["local_pooled_8bp"]==base["local_pooled_8bp"]})
    return {
        "evidence_class":EVIDENCE,"mode":"OFFLINE_ONLY", "assembly":"TOY_NOT_GRCH38",
        "source":{"filename":fixture_name,"sha256":sha256},"chrom":fixture["chrom"],
        "sequence_start_1based":fixture["sequence_start_1based"],
        "sequence_length_bp":len(fixture["sequence"]), "sequence":fixture["sequence"],
        "local_kmer":fixture["local_motif"],"distal_kmer":fixture["distal_motif"],
        "local_start_1based":fixture["local_motif_start_1based"],
        "distal_start_1based":fixture["distal_motif_start_1based"],
        "linear_distance_bp":fixture["distal_motif_start_1based"]-fixture["local_motif_start_1based"],
        "bin_size_bp":fixture["bin_size"],"bin_count":len(fixture["sequence"])//fixture["bin_size"],
        "cnn_analogy":"FIXED_EXACT_KMER_ONEHOT_SLIDING_FILTER_NO_TRAINING",
        "transformer_analogy":"HAND_SET_SINGLE_SCALAR_SOFTMAX_ATTENTION_NO_TRAINING",
        "decoder_analogy":"COARSE_CONTEXT_PLUS_LOCAL_SKIP_TO_1BP_TOY_OUTPUT",
        "cases":cases,"model_prediction":False,"official_alphagenome_architecture_reproduced":False,
        "dna_contact_map_measured":False,"clinical_interpretation":False,"real_api_calls_this_run":0,
        "disclaimer":"All sequences, attention weights and scores are artificial teaching proxies, not biology.",
    }


def csv_text(report: dict[str, Any]) -> str:
    buf=io.StringIO(newline="")
    out=csv.writer(buf,lineterminator="\n")
    out.writerow(["case","bin_index_0based","bin_start_1based","bin_end_1based_inclusive",
                  "cnn_local_full_hit","cnn_distal_full_hit","attention_from_local_bin","toy_fused_output"])
    for case in report["cases"]:
        a=case["local_context"]
        for i,(l,d,w) in enumerate(zip(a["local_pooled_8bp"],a["distal_pooled_8bp"],a["attention_from_local_bin"])):
            start=report["sequence_start_1based"]+i*report["bin_size_bp"]
            out.writerow([case["label"],i,start,start+report["bin_size_bp"]-1,l,d,f"{w:.8f}",f"{a['toy_fused_output']:.6f}"])
    return buf.getvalue()


def svg_text(report: dict[str, Any]) -> str:
    vals=[c["local_context"]["toy_fused_output"] for c in report["cases"]]
    d1=report["cases"][1]["delta_alt_minus_ref_toy_output"]
    d2=report["cases"][2]["delta_alt_minus_ref_toy_output"]
    colors={"bg":"#E5DEDA","paper":"#FAF8F6","ink":"#343642","muted":"#616577",
            "blue":"#639DB2","rose":"#C97E88","gold":"#B69B65","teal":"#669D95","line":"#D4C7C3"}
    def txt(x,y,msg,size=26,color=None,weight=400,anchor="start"):
        return f'<text x="{x}" y="{y}" font-family="Noto Sans CJK TC, sans-serif" font-size="{size}" fill="{color or colors["ink"]}" font-weight="{weight}" text-anchor="{anchor}">{html.escape(str(msg))}</text>'
    def card(x,y,w,h,fill=None,stroke=None,rx=25):
        return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill or colors["paper"]}" stroke="{stroke or colors["line"]}" stroke-width="2"/>'
    z=[]
    z.append('<svg xmlns="http://www.w3.org/2000/svg" width="1920" height="1080" viewBox="0 0 1920 1080" role="img" aria-label="Synthetic CNN and Transformer multiscale DNA learning diagram">')
    z.append(f'<rect width="1920" height="1080" fill="{colors["bg"]}"/>')
    z.append(txt(90,91,'Day 07  |  CNN × Transformer 多尺度學習',53,weight=700))
    z.append(txt(92,139,'從局部 3-mer 到跨區塊的資訊交換  ·  64 bp SYNTHETIC TOY, NOT ALPHAGENOME',23,colors["muted"]))
    # architecture row
    for x,w in [(85,520),(704,510),(1314,520)]:z.append(card(x,175,w,258))
    z.append(txt(125,225,'01  CNN-like 局部掃描',29,colors['blue'],700))
    z.append(txt(125,268,'固定 ACG / TGA 滑動比較核',22))
    z.append(txt(125,309,'每次看 3 個 DNA 字母',22))
    z.append(txt(125,355,'完全相同 → 局部命中 1',21,colors["muted"]))
    z.append(txt(125,397,'未訓練參數；不是正式 CNN',18,colors['rose']))
    z.append(txt(744,225,'02  下採樣（8 bp bins）',29,colors['gold'],700))
    z.append(txt(744,268,'64 bp → 8 個粗粒度 Token',22))
    z.append(txt(744,309,'局部命中：B1',22))
    z.append(txt(744,355,'遠端命中：B6',22))
    z.append(txt(744,397,'解析度改變，不等於訊號消失',18,colors['muted']))
    z.append(txt(1352,225,'03  Transformer-like 注意力',28,colors['rose'],700))
    z.append(txt(1352,268,'單一純量 Q / K / V + softmax',21))
    z.append(txt(1352,309,'讓 B1 向 B6 取得合成資訊',21))
    z.append(txt(1352,355,'再結合局部 skip 訊號',21))
    z.append(txt(1352,397,'手工權重；不是原模型',18,colors['muted']))
    z.append(f'<path d="M 620 303 L 687 303" stroke="{colors["blue"]}" stroke-width="7" stroke-linecap="round"/>')
    z.append(f'<path d="M 680 291 L 695 303 L 680 315" fill="none" stroke="{colors["blue"]}" stroke-width="5"/>')
    z.append(f'<path d="M 1230 303 L 1296 303" stroke="{colors["rose"]}" stroke-width="7" stroke-linecap="round"/>')
    z.append(f'<path d="M 1288 291 L 1305 303 L 1288 315" fill="none" stroke="{colors["rose"]}" stroke-width="5"/>')
    # sequence strip
    z.append(card(85,452,1749,276))
    z.append(txt(125,502,'TOY DNA：局部 ACG 與遠端 TGA 相隔 40 bp',29,weight=700))
    left=149; step=25.1; cy=595
    local_x=left+10*step+11; distant_x=left+50*step+11
    z.append(f'<path d="M {local_x:.1f} {cy-40} Q {(local_x+distant_x)/2:.1f} 465 {distant_x:.1f} {cy-40}" stroke="{colors["rose"]}" stroke-width="7" fill="none" stroke-dasharray="9 10"/>')
    z.append(txt((local_x+distant_x)/2,550,'Toy long-range link（非染色體接觸證據）',20,colors['rose'],700,"middle"))
    seq=report['sequence']
    for i,ch in enumerate(seq):
        is_local=10 <= i <=12
        is_distant=50<=i<=52
        fill=colors['blue'] if is_local else colors['rose'] if is_distant else '#F1EBE6'
        outline='#6892A5' if is_local else '#C28087' if is_distant else '#D7CFCC'
        x=left+i*step
        z.append(card(f'{x:.1f}',cy,21,37,fill,outline,5))
        z.append(txt(f'{x+10.5:.1f}',cy+26,ch,18,'#FFFFFF' if is_local or is_distant else colors['ink'],700,'middle'))
        if i in (0,10,50,63):z.append(txt(x+10,cy+65,str(report['sequence_start_1based']+i),17,colors['muted'],anchor='middle'))
    z.append(txt(127,700,'8 個 8 bp token',23,colors['muted']))
    for i in range(8):
        bx=400+i*149
        fill=('#DBEAF1' if i==1 else '#F2D8DC' if i==6 else '#EEE9E4')
        z.append(card(bx,661,115,42,fill,colors['line'],14))
        z.append(txt(bx+58,689,f'B{i}',23,colors['ink'],700,'middle'))
    # bottom bars
    z.append(txt(100,780,'不同人工改動對同一個 Toy Output 的影響',30,weight=700))
    cases=[('REF：局部 + 遠端',vals[0],None,colors['teal']),('遠端 T→G',vals[1],d1,colors['gold']),('局部 C→A',vals[2],d2,colors['rose'])]
    for i,(name,val,delta,clr) in enumerate(cases):
        x=85+i*586
        z.append(card(x,802,565,190))
        z.append(txt(x+30,850,name,25,weight=700))
        z.append(txt(x+530,851,f'{val:.4f}',29,clr,700,'end'))
        z.append(f'<rect x="{x+30}" y="874" width="505" height="30" rx="10" fill="#ECE8E5"/>')
        z.append(f'<rect x="{x+30}" y="874" width="{max(2,val*505):.1f}" height="30" rx="10" fill="{clr}"/>')
        z.append(txt(x+30,943,('REF 基準' if delta is None else f'ALT − REF = {delta:+.4f}'),21,colors['muted']))
        z.append(txt(x+30,974,'模型未訓練  ·  無生物學校準',16,colors['muted']))
    z.append(txt(89,1040,'RESEARCH EDUCATION ONLY  ·  NOT trained AlphaGenome  ·  NOT genomic causality  ·  API calls = 0',21,colors['muted']))
    z.append('</svg>')
    return '\n'.join(z)+'\n'


def html_text(report: dict[str, Any], svg: str) -> str:
    rows=[]
    for case in report['cases']:
        a=case['local_context']
        delta=case.get('delta_alt_minus_ref_toy_output')
        rows.append('<tr><td>'+html.escape(case['label'])+'</td><td>'+html.escape(case.get('variant_id') or 'REF')+'</td>'
                    +f'<td>{len(a["local_kernel_full_hits_1based"])}</td><td>{len(a["distal_kernel_full_hits_1based"])}</td>'
                    +f'<td>{a["toy_fused_output"]:.6f}</td><td>{"—" if delta is None else f"{delta:+.6f}"}</td></tr>')
    return ('''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Day07 CNN Transformer | Synthetic Lab</title><style>body{margin:0;background:#E5DEDA;color:#323640;font-family:system-ui,'Noto Sans CJK TC',sans-serif}.wrap{max-width:1160px;margin:auto;padding:30px}h1{font-size:2rem}figure{margin:10px 0}svg{width:100%;height:auto;border-radius:16px}table{border-collapse:collapse;width:100%;background:#fffdf9}th,td{padding:14px;border-bottom:1px solid #ddd;text-align:left}th{background:#e9e7e7}.box{background:#fffaf7;padding:22px;border-radius:16px}code{word-break:break-all}small{color:#615f63}</style></head><body><main class="wrap"><h1>Day07｜CNN 與 Transformer 多尺度教學代理</h1><p><strong>SYNTHETIC / OFFLINE ONLY</strong>：固定 k-mer 掃描與手工純量 attention；並非 AlphaGenome 原模型。</p><figure>'''
            +svg+'''</figure><section class="box"><h2>比較三組人工序列</h2><table><thead><tr><th>情境</th><th>Variant</th><th>CNN local hits</th><th>Distal hits</th><th>Toy output</th><th>ALT−REF</th></tr></thead><tbody>'''
            +''.join(rows)+'''</tbody></table><details><summary>為何遠端改動會改變本地分數？</summary><p>本例手動設計 scalar attention，把遠端 motif 命中轉成數學上下文。這不是已發現的 enhancer–promoter 作用，更不是 3D contact map。</p></details></section><p><small>資料類別：SYNTHETIC_EDUCATIONAL_DEMO。無患者資料，無 API，沒有訓練權重。SHA-256 可檢查本地檔案一致性，不證明來源真實。</small></p></main></body></html>\n''')


def render_expected(report: dict[str, Any]) -> dict[str, bytes]:
    svg=svg_text(report)
    return {
        OUTPUT_FILES[0]:(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+'\n').encode('utf-8'),
        OUTPUT_FILES[1]:csv_text(report).encode('utf-8'),
        OUTPUT_FILES[2]:svg.encode('utf-8'),
        OUTPUT_FILES[3]:html_text(report,svg).encode('utf-8'),
    }


def manifest_text(report: dict[str, Any], files: dict[str, bytes]) -> str:
    manifest={"evidence_class":EVIDENCE,"assembly":"TOY_NOT_GRCH38",
              "source":report['source'],"hash_algorithm":"SHA-256",
              "output_sha256":{name:hashlib.sha256(data).hexdigest() for name,data in sorted(files.items())},
              "signed_provenance":False,"official_model_output":False,"real_api_calls_this_run":0}
    return json.dumps(manifest,ensure_ascii=False,sort_keys=True,indent=2)+'\n'


def run(fixture_path: str | Path, output_dir: str | Path) -> dict[str, Any]:
    fixture, sha=load_fixture(fixture_path)
    report=build_report(fixture,sha,Path(fixture_path).name)
    expected=render_expected(report)
    directory=Path(output_dir)
    directory.mkdir(parents=True, exist_ok=False)  # Never overwrite existing work.
    for name,payload in expected.items():
        (directory/name).write_bytes(payload)
    (directory/MANIFEST_NAME).write_text(manifest_text(report,expected),encoding='utf-8')
    cases=report['cases']
    return {"status":"PASS","mode":"OFFLINE_ONLY","case_count":len(cases),
            "local_motif_start_1based":report['local_start_1based'],
            "distal_motif_start_1based":report['distal_start_1based'],
            "distance_bp":report['linear_distance_bp'],
            "ref_toy_output":cases[0]['local_context']['toy_fused_output'],
            "distal_variant_toy_output":cases[1]['local_context']['toy_fused_output'],
            "local_variant_toy_output":cases[2]['local_context']['toy_fused_output'],
            "artifact_count":len(expected),"real_api_calls_this_run":0}


def verify(fixture_path: str | Path, output_dir: str | Path) -> dict[str, Any]:
    fixture,sha=load_fixture(fixture_path)
    expected_report=build_report(fixture,sha,Path(fixture_path).name)
    expected=render_expected(expected_report)
    folder=Path(output_dir)
    actual_manifest=json.loads((folder/MANIFEST_NAME).read_text(encoding='utf-8'))
    expected_manifest=json.loads(manifest_text(expected_report,expected))
    if actual_manifest != expected_manifest:
        raise ValueError("Manifest or source is inconsistent with expected toy report")
    for name, data in expected.items():
        stored=(folder/name).read_bytes()
        if hashlib.sha256(stored).hexdigest()!=actual_manifest['output_sha256'][name]:
            raise ValueError("SHA mismatch: "+name)
        if stored!=data:
            raise ValueError("Content mismatch: "+name)
    return {"status":"PASS","checksums":"MATCH","output_artifacts_verified":len(expected),
            "source_sha256_verified":True,"real_api_calls_this_run":0}


def main(argv: list[str] | None=None) -> int:
    parser=argparse.ArgumentParser(description="Offline synthetic CNN/Transformer multiscale lesson")
    parser.add_argument('command',choices=('run','verify'))
    parser.add_argument('--mode',choices=('offline',),default='offline')
    parser.add_argument('--fixture',default='day07_synthetic_multiscale.json')
    parser.add_argument('--output-dir',default='day07_output')
    a=parser.parse_args(argv)
    result=(run(a.fixture,a.output_dir) if a.command=='run' else verify(a.fixture,a.output_dir))
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
