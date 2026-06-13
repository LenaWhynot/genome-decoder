"""Импутированные данные Michigan → нормализованный файл для пайплайна (шаги 4/5/13).

Из 23 chr*.dose.vcf.gz берёт колонку образца (1-я из 5 дублей), фильтрует по качеству
импутации (R²≥0.3 ИЛИ TYPED=реально типизировано), оставляет сайты, где она несёт ALT,
и пишет в формате step1 (rsid,chrom,pos,a1,a2,zygosity,ambiguous) → out/01_normalized.tsv.

Так шаг 4 (ClinVar по позиции) и далее увидят ~миллионы позиций вместо 592К.
ID импутированных = chr:pos (rsID нет) — для ClinVar-матча по позиции это неважно.
Фарма/трейты/нутри (по rsID) остаются на типизированных (01_normalized_typed.tsv).
"""
import gzip
from pathlib import Path
import config

IMP = config.BASE / ".cache" / "imputed"
R2_MIN = 0.3
AMB = config.AMBIGUOUS_PAIRS


def parse_info(info):
    d = {}
    for f in info.split(";"):
        if "=" in f:
            k, v = f.split("=", 1)
            d[k] = v
        else:
            d[f] = True
    return d


def run():
    chroms = [str(c) for c in range(1, 23)] + ["X"]
    n_total = n_pass = n_carry = 0
    out = open(config.NORMALIZED, "w", encoding="utf-8")
    out.write("rsid\tchrom\tpos\ta1\ta2\tzygosity\tambiguous\tr2\ttyped\n")

    for ch in chroms:
        fp = IMP / f"chr{ch}.dose.vcf.gz"
        if not fp.exists():
            print(f"[imp] нет {fp.name}, пропуск")
            continue
        c_pass = c_carry = 0
        with gzip.open(fp, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                n_total += 1
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10:
                    continue
                chrom, pos, vid, ref, alt, _q, _flt, info, fmt, s1 = cols[:10]
                inf = parse_info(info)
                typed = "TYPED" in inf
                try:
                    r2 = float(inf.get("R2", "0"))
                except ValueError:
                    r2 = 0.0
                if not (typed or r2 >= R2_MIN):
                    continue
                c_pass += 1
                gt = s1.split(":")[0].replace("|", "/")
                if "1" not in gt:           # не несёт ALT
                    continue
                # только биаллельные SNV (ref/alt одной буквы) — для матча с ClinVar/частотами
                if len(ref) != 1 or len(alt) != 1:
                    continue
                a = gt.split("/")
                a1 = ref if a[0] == "0" else alt
                a2 = ref if a[1] == "0" else alt
                zyg = "hom" if a1 == a2 else "het"
                amb = 1 if frozenset({a1, a2}) in AMB else 0
                r2s = f"{r2:.3f}" if not typed else "1.0"
                out.write(f"{vid}\t{chrom}\t{pos}\t{a1}\t{a2}\t{zyg}\t{amb}\t{r2s}\t{int(typed)}\n")
                c_carry += 1
        n_pass += c_pass
        n_carry += c_carry
        print(f"[imp] chr{ch}: прошло R²/typed {c_pass:,}, несёт ALT (SNV) {c_carry:,}")

    out.close()
    print(f"\n[imp] всего вариантов в dose: {n_total:,}")
    print(f"[imp] прошло качество (R²≥{R2_MIN} или typed): {n_pass:,}")
    print(f"[imp] несёт ALT, биаллельный SNV → в нормализованный: {n_carry:,}")
    print(f"[imp] → {config.NORMALIZED}")


if __name__ == "__main__":
    run()
