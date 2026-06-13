"""ШАГ 5b — gnomAD-частоты напрямую (без OpenCRAVAT/аккаунта/80 ГБ).

Запрашивает популяционные частоты gnomAD r2.1 (GRCh37, совпадает с данными образца) ТОЛЬКО
для значимых вариантов из шага 5 — через публичный GraphQL-API, батчами. Это «настоящий
частотный фильтр здравого смысла» Григория: AF≥5% → доброкачественно по частоте (ACMG BA1),
подтверждает остаточные гомозиготные артефакты как частые/неопасные.
"""
import json
import time
import urllib.request
import config

API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r2_1"          # GRCh37
BATCH = 20
OUT = config.OUT / "05b_gnomad.tsv"


def load_significant():
    rows = []
    with open(config.SIGNIFICANT, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(h, line.rstrip("\n").split("\t"))))
    return rows


def _popmax(node):
    """Макс AF по global + субпопуляциям (BA1-консервативно). af поп = ac/an."""
    if not node:
        return []
    afs = []
    if node.get("af") is not None:
        afs.append(node["af"])
    for p in (node.get("populations") or []):
        an = p.get("an") or 0
        ac = p.get("ac")
        if ac is not None and an >= 2000:        # достаточный размер выборки
            afs.append(ac / an)
    return afs


def query_batch(variant_ids):
    """variant_ids → {vid: max_af or None}. Берём popmax (не только global), чек errors."""
    parts = []
    for i, vid in enumerate(variant_ids):
        parts.append(
            f'v{i}: variant(variantId: "{vid}", dataset: {DATASET}) '
            f'{{ genome {{ af populations {{ id ac an }} }} '
            f'exome {{ af populations {{ id ac an }} }} }}'
        )
    query = "query {\n" + "\n".join(parts) + "\n}"
    body = json.dumps({"query": query}).encode()
    req = urllib.request.Request(API, data=body,
                                 headers={"Content-Type": "application/json"})
    data = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read())
            break
        except Exception as e:
            if attempt == 3:
                print(f"[step5b]   батч ПРОПУЩЕН (сеть): {e}")
                return {vid: "ERR" for vid in variant_ids}  # не путать с «редкий»
            time.sleep(3 * (attempt + 1))
    out = {}
    d = (data or {}).get("data") or {}
    # ошибки GraphQL по части алиасов — НЕ выдавать за «редкий»
    err_paths = set()
    for e in (data or {}).get("errors") or []:
        for seg in (e.get("path") or []):
            if isinstance(seg, str) and seg.startswith("v"):
                err_paths.add(seg)
    for i, vid in enumerate(variant_ids):
        key = f"v{i}"
        v = d.get(key)
        if v is None and key in err_paths:
            out[vid] = "ERR"            # API-ошибка, не «нет в gnomAD»
            continue
        afs = _popmax(v.get("genome") if v else None) + _popmax(v.get("exome") if v else None)
        out[vid] = max(afs) if afs else None
    return out


def run():
    rows = load_significant()
    # уникальные variantId (chrom-pos-ref-alt, без 'chr')
    ids = []
    seen = set()
    for r in rows:
        vid = f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}"
        if vid not in seen and len(r['ref']) == 1 and len(r['alt']) == 1:
            seen.add(vid)
            ids.append(vid)
    print(f"[step5b] значимых вариантов на запрос в gnomAD: {len(ids)}")

    af_map = {}
    for i in range(0, len(ids), BATCH):
        batch = ids[i:i + BATCH]
        af_map.update(query_batch(batch))
        print(f"[step5b]   {min(i + BATCH, len(ids))}/{len(ids)}")
        time.sleep(1)

    with open(OUT, "w", encoding="utf-8") as out:
        out.write("category\tgene\trsid\tzygosity\tclnsig\tstars\tgnomad_af\tverdict\n")
        for r in rows:
            vid = f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}"
            af = af_map.get(vid)
            if af == "ERR":
                verdict = "gnomAD: ошибка API (не оценено)"
                afs = ""
            elif af is None:
                verdict = "не в gnomAD (редкий/н-д)"
                afs = ""
            else:
                afs = f"{af:.5f}"
                if af >= 0.05:
                    verdict = f"ЧАСТЫЙ {af:.1%} (popmax) → доброкач. по частоте (BA1)"
                elif af >= 0.01:
                    verdict = f"низкопенетр. {af:.1%} (popmax)"
                else:
                    verdict = f"редкий {af:.3%} (popmax)"
            out.write(f"{r['category']}\t{r['gene']}\t{r['rsid']}\t{r['zygosity']}\t"
                      f"{r['clnsig']}\t{r['stars']}\t{afs}\t{verdict}\n")

    # сводка: что gnomAD подтверждает как частое (артефакты) vs реально редкое
    plp = [r for r in rows if r["category"] == "P/LP"]
    common = rare = unknown = err = 0
    for r in plp:
        af = af_map.get(f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}")
        if af == "ERR":
            err += 1
        elif af is None:
            unknown += 1
        elif af >= 0.05:
            common += 1
        else:
            rare += 1
    if err:
        print(f"[step5b]   ⚠ API-ошибок (не оценено): {err}")
    print(f"[step5b] P/LP всего: {len(plp)}")
    print(f"[step5b]   gnomAD ЧАСТЫЙ ≥5% (→ артефакт/доброкач.): {common}")
    print(f"[step5b]   реально редкий <5%: {rare}")
    print(f"[step5b]   нет в gnomAD (ультраредкий или н-д): {unknown}")
    print(f"[step5b] → {OUT}")


if __name__ == "__main__":
    run()
