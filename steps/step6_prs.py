"""ШАГ 6 — Полигенные шкалы риска (PRS) из PGS Catalog.

Для каждой болезни: scoring-файл (вариант→эффект-аллель→вес) из PGS Catalog (GRCh37,
гармонизированный), берёт дозу эффект-аллеля образца из импутированных dose-файлов,
score = Σ weight·dosage. Перцентиль — аналитически через частоты (HWE): где образец в
распределении популяции.

Оговорки: (1) допущение независимости вариантов (без LD) → дисперсия завышена,
перцентиль КОНСЕРВАТИВЕН (крайности занижены); (2) частоты из шкалы (часто европейские);
(3) импутированные дозы. НЕ диагноз — относительный риск, обсуждать с врачом.
"""
import gzip
import math
import urllib.request
from pathlib import Path
import config

PGS = {
    "PGS000004": "Рак груди (Mavaddat 313)",
    "PGS000349": "Ишемическая болезнь сердца",
    "PGS000043": "Венозная тромбоэмболия (тромбозы)",
    "PGS000037": "Астма",
}
IMP = config.BASE / ".cache" / "imputed"
PGSDIR = config.CACHE / "pgs"
PGSDIR.mkdir(exist_ok=True)
R2_MIN = 0.3
OUT = config.OUT / "06_prs.tsv"


def phi(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def download(pgs):
    fp = PGSDIR / f"{pgs}.txt.gz"
    if fp.exists() and fp.stat().st_size > 1000:
        return fp
    url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs}/ScoringFiles/Harmonized/{pgs}_hmPOS_GRCh37.txt.gz"
    print(f"[step6] качаю {pgs} …")
    urllib.request.urlretrieve(url, fp)
    return fp


def parse_pgs(fp):
    """{(chrom,pos): (effect, other, weight, af)} — только с частотой."""
    variants = {}
    with gzip.open(fp, "rt", encoding="utf-8") as f:
        cols = None
        for line in f:
            if line.startswith("#"):
                continue
            if cols is None:
                cols = line.rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(cols)}
                continue
            p = line.rstrip("\n").split("\t")
            try:
                chrom = p[idx["hm_chr"]]
                pos = p[idx["hm_pos"]]
                eff = p[idx["effect_allele"]].upper()
                oth = p[idx.get("other_allele", idx.get("hm_inferOtherAllele"))].upper()
                w = float(p[idx["effect_weight"]])
                afi = idx.get("allelefrequency_effect")
                af = float(p[afi]) if afi is not None and p[afi] not in ("", "NA") else None
            except (KeyError, ValueError, IndexError):
                continue
            if not chrom or not pos or af is None or not (0 < af < 1):
                continue
            # список на позицию — multiallelic не перетирается
            variants.setdefault((chrom, pos), []).append((eff, oth, w, af))
    return variants


def build_dosage(union):
    """{(chrom,pos): (ref,alt,ds)} из импутированных dose-файлов для нужных позиций."""
    chroms = sorted({c for c, _ in union}, key=lambda x: (x.isdigit(), x))
    dos = {}
    for ch in [str(c) for c in range(1, 23)] + ["X"]:
        fp = IMP / f"chr{ch}.dose.vcf.gz"
        if not fp.exists():
            continue
        want = {pos for (c, pos) in union if c == ch}
        if not want:
            continue
        with gzip.open(fp, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10 or cols[1] not in want:
                    continue
                chrom, pos, _id, ref, alt, _q, _flt, info, fmt, s1 = cols[:10]
                # R²-фильтр НЕ применяем: трюк с 5 дублями занизил R², но дозы валидны
                ff = fmt.split(":")
                if "DS" not in ff:
                    continue
                try:
                    ds = float(s1.split(":")[ff.index("DS")])
                except (ValueError, IndexError):
                    continue
                dos.setdefault((chrom, pos), []).append((ref.upper(), alt.upper(), ds))
    return dos


def run():
    parsed = {}
    union = set()
    for pgs in PGS:
        fp = download(pgs)
        v = parse_pgs(fp)
        parsed[pgs] = v
        union |= set(v.keys())
        print(f"[step6] {pgs}: вариантов с частотой {sum(len(x) for x in v.values()):,}")

    print(f"[step6] уникальных позиций всего: {len(union):,} — собираю дозы из импутации…")
    dos = build_dosage(union)
    print(f"[step6] нашлось доз в импутации: {len(dos):,}")

    rows = []
    for pgs, variants in parsed.items():
        score = mu = var = 0.0
        matched = mism = 0
        total = sum(len(v) for v in variants.values())   # все варианты (с multiallelic)
        for (chrom, pos), vlist in variants.items():
            dlist = dos.get((chrom, pos))
            if not dlist:
                continue
            for (eff, oth, w, af) in vlist:
                # среди доз на этой позиции найти запись с тем же набором аллелей
                eff_dos = None
                for (ref, alt, ds) in dlist:
                    if {eff, oth} == {ref, alt}:
                        eff_dos = ds if eff == alt else 2 - ds
                        break
                if eff_dos is None:
                    mism += 1
                    continue
                matched += 1
                score += w * eff_dos
                mu += w * 2 * af
                var += (w ** 2) * 2 * af * (1 - af)
        if not total or var <= 0:
            rows.append((pgs, PGS[pgs], matched, total, 0.0, score, None, None,
                         "н/д — scoring-файл без частот популяции"))
            continue
        z = (score - mu) / math.sqrt(var)
        pct = 100 * phi(z)
        band = ("высокий (>95%)" if pct >= 95 else
                "повышенный (80–95%)" if pct >= 80 else
                "ниже среднего (<20%)" if pct < 20 else "средний (20–80%)")
        cover = matched / total * 100
        rows.append((pgs, PGS[pgs], matched, total, cover, score, z, pct, band))

    with open(OUT, "w", encoding="utf-8") as out:
        out.write("pgs\ttrait\tmatched\ttotal\tcoverage_pct\traw_score\tz\tpercentile\tband\n")
        for r in rows:
            pcts = f"{r[7]:.0f}" if r[7] is not None else ""
            zs = f"{r[6]:.2f}" if r[6] is not None else ""
            out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.0f}\t{r[5]:.4f}\t{zs}\t{pcts}\t{r[8]}\n")

    print("\n[step6] РЕЗУЛЬТАТ (перцентиль = где образец в популяции):")
    for r in rows:
        if r[7] is not None:
            print(f"[step6]   {r[1]}: ~{r[7]:.0f}-й перцентиль · {r[8]} "
                  f"(покрытие {r[4]:.0f}%, {r[2]}/{r[3]})")
        else:
            print(f"[step6]   {r[1]}: {r[8]} (нужна референсная панель частот)")
    print(f"[step6] → {OUT}")
    return rows


if __name__ == "__main__":
    run()
