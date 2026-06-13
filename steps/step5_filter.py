"""ШАГ 5 — Фильтр значимости.

Аналог фильтра здравого смысла Григория (у него — частота gnomAD; полная база gnomAD
слишком тяжёлая для локали, поэтому в роли фильтра доверия используем уровень ClinVar —
«золотые звёзды» review status: сколько лабораторий подтвердили).

Оставляем то, где образец несёт оценённый аллель И оценка клинически значима.
Категории: P/LP (патогенное), risk (фактор риска/ассоциация), pgx (ответ на лекарства),
conflict (конфликтующие оценки — на ручную проверку).
"""
import config

SEVERITY = {"P/LP": 0, "risk": 1, "pgx": 2, "conflict": 3}


def categorize(clnsig: str):
    s = clnsig.lower()
    if "conflict" in s:
        return "conflict"
    if "drug response" in s:
        return "pgx"
    if "pathogenic" in s:           # ловит Pathogenic, Likely pathogenic, Pathogenic/Likely
        return "P/LP"
    if "risk factor" in s or "association" in s or "affects" in s:
        return "risk"
    return None  # Benign / Uncertain / прочее — не значимо


def freq_flag(max_af: str):
    """Частотный фильтр здравого смысла (шаг 5 Григория) по зашитым в ClinVar частотам.
    Возвращает (drop?, flag). AF≥5% = доброкач. по частоте (ACMG BA1)."""
    if not max_af:
        return False, "частота н/д"
    af = float(max_af)
    if af >= 0.05:
        return True, f"частый {af:.1%} (отсев)"
    if af >= 0.01:
        return False, f"низкопенетр. {af:.1%}"
    return False, f"редкий {af:.2%}"


def run():
    rows = []
    benign = uncertain = freq_dropped = 0
    with open(config.CLINVAR_HITS, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            vals = line.rstrip("\n").split("\t")
            rec = dict(zip(header, vals))
            cat = categorize(rec["clnsig"])
            if cat is None:
                s = rec["clnsig"].lower()
                if "benign" in s:
                    benign += 1
                elif "uncertain" in s:
                    uncertain += 1
                continue
            drop, flag = freq_flag(rec.get("max_af", ""))
            if drop:
                freq_dropped += 1
                continue
            rec["freq_flag"] = flag
            rec["category"] = cat
            rec["stars"] = int(rec["stars"])
            rec["alt_copies"] = int(rec["alt_copies"])
            rows.append(rec)

    # сортировка: сначала категория (P/LP важнее), потом доверие (звёзды), потом гомозиготы
    rows.sort(key=lambda r: (SEVERITY[r["category"]], -r["stars"], -r["alt_copies"]))

    cols = ["category", "stars", "freq_flag", "max_af", "clnsig", "gene", "rsid", "zygosity",
            "my_a1", "my_a2", "ref", "alt", "ambiguous", "condition", "chrom", "pos", "clnvc",
            "r2", "typed"]
    with open(config.SIGNIFICANT, "w", encoding="utf-8") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    from collections import Counter
    by_cat = Counter(r["category"] for r in rows)
    hi = sum(1 for r in rows if r["stars"] >= 2)
    hi_rare = sum(1 for r in rows if r["stars"] >= 2 and r["category"] == "P/LP"
                  and not r["freq_flag"].startswith("низкопенетр"))
    print(f"[step5] значимых находок (несёт аллель): {len(rows):,}")
    print(f"[step5]   по категориям: {dict(by_cat)}")
    print(f"[step5]   с доверием ≥2★: {hi} (из них P/LP редких ≥2★: {hi_rare})")
    print(f"[step5]   отсеяно по частоте ≥5% (BA1): {freq_dropped:,}")
    print(f"[step5]   отсеяно Benign: {benign:,} · Uncertain: {uncertain:,}")
    print(f"[step5] → {config.SIGNIFICANT}")
    return rows


if __name__ == "__main__":
    run()
