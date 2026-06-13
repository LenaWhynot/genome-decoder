"""ШАГ 8 (AI-скор) — AlphaMissense (DeepMind), account-free.

Без OpenCRAVAT/аккаунта: качаем готовые предсказания AlphaMissense hg19 (Zenodo, 622МБ,
лежат в .cache/) и аннотируем миссенс-варианты, которые образец РЕАЛЬНО несёт. Это «4-я
нейросеть голосует опасен/нет» из слайда — особенно ценно для вариантов, которых нет
в ClinVar (наука их ещё не классифицировала).

Вход — Genotek VCF (hg19, есть REF/ALT/GT, chr-префикс). Совпадает по цепи с AlphaMissense.
"""
import gzip
from pathlib import Path

BASE = Path(__file__).resolve().parent
IMP = BASE / ".cache" / "imputed"
AM_GZ = BASE / ".cache" / "AlphaMissense_hg19.tsv.gz"
OUT = BASE / "out" / "08_alphamissense.tsv"
OUT.parent.mkdir(exist_ok=True)
R2_MIN = 0.3


def load_carried():
    """{(chrom,pos): {(ref,alt): 'het'/'hom'}} — варианты на ИМПУТИРОВАННЫХ данных, где
    образец несёт ALT (R²≥0.3 или typed). chrom приводим к 'chr1' под формат AlphaMissense.
    Раньше брали типизированный VCF — правка по ревью (ClinVar на импутированных → AM тоже)."""
    carried = {}
    for ch in [str(c) for c in range(1, 23)] + ["X"]:
        fp = IMP / f"chr{ch}.dose.vcf.gz"
        if not fp.exists():
            continue
        with gzip.open(fp, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 10:
                    continue
                chrom, pos, _id, ref, alt, _q, _flt, info, fmt, s1 = c[:10]
                r2 = 0.0
                for fld in info.split(";"):
                    if fld.startswith("R2="):
                        try:
                            r2 = float(fld[3:])
                        except ValueError:
                            pass
                if "TYPED" not in info and r2 < R2_MIN:
                    continue
                # multiallelic: ALT может быть "A,C", GT — 0/2, 1/2, 1/1…
                alts = alt.split(",")
                gt = s1.split(":")[0].replace("|", "/")
                idxs = [a for a in gt.split("/") if a.isdigit()]
                carried_idx = {int(a) for a in idxs if int(a) > 0}   # какие ALT-аллели несёт
                for ai in carried_idx:
                    if ai - 1 >= len(alts):
                        continue
                    al = alts[ai - 1]
                    if len(ref) != 1 or len(al) != 1:   # только SNV
                        continue
                    zyg = "hom" if idxs.count(str(ai)) == 2 else "het"
                    carried.setdefault(("chr" + chrom, pos), {})[(ref, al)] = zyg
    return carried


def run():
    carried = load_carried()
    print(f"[step8] позиций, где образец несёт ALT: {len(carried):,}")

    n_scanned = n_hit = 0
    rows = []
    with gzip.open(AM_GZ, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_scanned += 1
            chrom, pos, ref, alt, genome, uni, tx, pvar, score, cls = line.rstrip("\n").split("\t")
            key = (chrom, pos)
            hit = carried.get(key)
            if not hit:
                continue
            zyg = hit.get((ref, alt))
            if zyg is None:
                continue
            n_hit += 1
            rows.append((chrom, pos, ref, alt, pvar, float(score), cls, zyg, uni, tx))

    # сортировка: сначала самые «патогенные» по скору
    rows.sort(key=lambda r: -r[5])
    with open(OUT, "w", encoding="utf-8") as out:
        out.write("chrom\tpos\tref\talt\tprotein_variant\tam_pathogenicity\tam_class\tzygosity\tuniprot\ttranscript\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    from collections import Counter
    by_cls = Counter(r[6] for r in rows)
    patho = [r for r in rows if r[6] == "pathogenic"]
    print(f"[step8] строк AlphaMissense просмотрено: {n_scanned:,}")
    print(f"[step8] миссенс-вариантов образца со скором: {n_hit:,}")
    print(f"[step8]   по классам: {dict(by_cls)}")
    print(f"[step8]   likely_pathogenic (топ по скору):")
    for r in patho[:12]:
        print(f"[step8]     {r[8]} {r[4]} score={r[5]:.3f} {r[7]} (chr{r[0].replace('chr','')}:{r[1]})")
    print(f"[step8] → {OUT}")
    return rows


if __name__ == "__main__":
    run()
