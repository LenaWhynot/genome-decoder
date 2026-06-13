"""ШАГ 8 (доп) — CADD + SpliceAI для значимых вариантов (API, GRCh37).

Добивает «4 нейросети»: AlphaMissense уже есть (step8), тут добавляем CADD (общая
вредность, PHRED) и SpliceAI (нарушение сплайсинга) — точечно по значимым вариантам,
через публичные API (без 80 ГБ локально). SpliceAI-API нестабилен → best-effort.
"""
import json
import time
import urllib.request
import config

CADD = "https://cadd.gs.washington.edu/api/v1.0/GRCh37-v1.6/"
SPLICEAI = "https://spliceailookup-api.broadinstitute.org/spliceai/"
OUT = config.OUT / "08b_cadd_spliceai.tsv"


def load_significant():
    rows = []
    with open(config.SIGNIFICANT, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            if len(r["ref"]) == 1 and len(r["alt"]) == 1:
                rows.append(r)
    return rows


def cadd_phred(chrom, pos, ref, alt):
    try:
        with urllib.request.urlopen(f"{CADD}{chrom}:{pos}", timeout=20) as r:
            data = json.loads(r.read())
        for rec in data:
            if rec.get("Alt") == alt and rec.get("Ref") == ref:   # сверяем и REF (правка ревью)
                return rec.get("PHRED", "")
    except Exception:
        return ""
    return ""


# SpliceAI lookup API стабильно отдаёт пусто (origin issue) → отключён. Сплайс-эффект
# по значимым уже частично закрыт VEP (step3: splice_acceptor/donor консеквенсы).
SPLICEAI_ENABLED = False


def spliceai_max(chrom, pos, ref, alt):
    if not SPLICEAI_ENABLED:
        return ""
    try:
        url = f"{SPLICEAI}?hg=37&distance=50&mask=0&variant={chrom}-{pos}-{ref}-{alt}"
        with urllib.request.urlopen(url, timeout=6) as r:
            data = json.loads(r.read())
        scores = data.get("scores") or []
        if scores and isinstance(scores[0], dict):
            ds = [scores[0].get(k, 0) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")]
            ds = [float(x) for x in ds if x not in ("", None)]
            return max(ds) if ds else ""
    except Exception:
        return ""
    return ""


def run():
    rows = load_significant()
    print(f"[step8b] значимых SNV: {len(rows)} — запрашиваю CADD (+SpliceAI best-effort)…")
    out_rows = []
    spliceai_ok = 0
    for i, r in enumerate(rows):
        phred = cadd_phred(r["chrom"], r["pos"], r["ref"], r["alt"])
        sai = spliceai_max(r["chrom"], r["pos"], r["ref"], r["alt"])
        if sai != "":
            spliceai_ok += 1
        out_rows.append((r["gene"], r["rsid"], r["category"], phred, sai, r["clnsig"]))
        if (i + 1) % 20 == 0:
            print(f"[step8b]   {i + 1}/{len(rows)}")
        time.sleep(0.3)

    with open(OUT, "w", encoding="utf-8") as out:
        out.write("gene\trsid\tcategory\tcadd_phred\tspliceai_max\tclnsig\n")
        for r in out_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    high_cadd = [r for r in out_rows if r[3] and _f(r[3]) >= 20]
    print(f"[step8b] CADD PHRED≥20 (топ-1% вредных): {len(high_cadd)}")
    for r in sorted(high_cadd, key=lambda x: -_f(x[3]))[:10]:
        print(f"[step8b]   {r[0]} {r[1]} CADD={r[3]} {('SpliceAI=' + str(r[4])) if r[4] else ''} ({r[2]})")
    print(f"[step8b] SpliceAI получено для {spliceai_ok}/{len(rows)} (API нестабилен)")
    print(f"[step8b] → {OUT}")
    return out_rows


def _f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0.0


if __name__ == "__main__":
    run()
