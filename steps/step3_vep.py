"""ШАГ 3 — Аннотация эффекта (Ensembl VEP REST API, GRCh37).

Без локального VEP/Docker: для значимых вариантов (шаг 5) спрашиваем Ensembl, ЧТО вариант
делает с белком — most_severe_consequence (missense / stop_gained / splice / regulatory…)
+ impact (HIGH/MODERATE/LOW) + ген. Это и есть «разметка эффекта» из транскрипта.
Только по значимым (не весь геном) — точечно, без 80 ГБ.
"""
import json
import time
import urllib.request
import config

VEP = "https://grch37.rest.ensembl.org/vep/human/region"
OUT = config.OUT / "03_vep.tsv"
BATCH = 150


def load_significant():
    rows = []
    with open(config.SIGNIFICANT, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(h, line.rstrip("\n").split("\t"))))
    return rows


def query_vep(variants):
    """variants: ['6 26093141 . G A . . .', ...] → {input_str: (consequence, impact, gene)}."""
    body = json.dumps({"variants": variants}).encode()
    req = urllib.request.Request(VEP, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                data = json.loads(r.read())
            break
        except Exception as e:
            if attempt == 3:
                print(f"[step3] батч пропущен: {e}")
                return {}
            time.sleep(3 * (attempt + 1))
    out = {}
    for v in data:
        inp = v.get("input", "")
        cons = v.get("most_severe_consequence", "")
        impact = ""
        gene = ""
        tcs = v.get("transcript_consequences") or []
        # берём транскрипт с самым сильным impact
        rank = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
        best = None
        for tc in tcs:
            if cons in (tc.get("consequence_terms") or []):
                if best is None or rank.get(tc.get("impact", ""), 0) > rank.get(best.get("impact", ""), 0):
                    best = tc
        if best:
            impact = best.get("impact", "")
            gene = best.get("gene_symbol", "")
        out[inp] = (cons, impact, gene)
    return out


def run():
    rows = load_significant()
    # уникальные SNV-варианты
    seen = {}
    for r in rows:
        if len(r["ref"]) == 1 and len(r["alt"]) == 1:
            vstr = f"{r['chrom']} {r['pos']} . {r['ref']} {r['alt']} . . ."
            seen[vstr] = r
    variants = list(seen.keys())
    print(f"[step3] значимых вариантов на аннотацию VEP: {len(variants)}")

    ann = {}
    for i in range(0, len(variants), BATCH):
        ann.update(query_vep(variants[i:i + BATCH]))
        print(f"[step3]   {min(i + BATCH, len(variants))}/{len(variants)}")
        time.sleep(1)

    from collections import Counter
    cons_count = Counter()
    with open(OUT, "w", encoding="utf-8") as out:
        out.write("gene\trsid\tcategory\tconsequence\timpact\tclnsig\n")
        for vstr, r in seen.items():
            cons, impact, gene = ann.get(vstr, ("", "", ""))
            cons_count[cons or "?"] += 1
            out.write(f"{gene or r['gene']}\t{r['rsid']}\t{r['category']}\t{cons}\t{impact}\t{r['clnsig']}\n")

    high = [v for v in seen.values()]
    print(f"[step3] эффекты (most_severe): {dict(cons_count)}")
    print(f"[step3] → {OUT}")
    return cons_count


if __name__ == "__main__":
    run()
