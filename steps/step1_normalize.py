"""ШАГ 1 — Нормализация. Авто-детект формата входа → единый TSV.

Поддержка:
  • 23andMe / AncestryDNA raw .txt  (TAB: rsid, chromosome, position, genotype)
  • VCF (.vcf / .vcf.gz)            (CHROM POS ID REF ALT ... GT) — raw или импутированный

Выход: NORMALIZED_TYPED (rsid, chrom без 'chr', pos, a1, a2, zygosity, ambiguous).
Отсев no-call; флаг палиндромных A/T,C/G (pre-QC). Сборка build 37 (как у raw-чипов).
"""
import gzip
import config


def _open(p):
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if str(p).endswith(".gz") \
        else open(p, encoding="utf-8", errors="replace")


def detect_format(path):
    with _open(path) as f:
        for line in f:
            if line.startswith("##fileformat=VCF") or line.startswith("#CHROM\t"):
                return "vcf"
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 4 and len(parts[3]) <= 2:
                return "23andme"
            if len(parts) >= 8:
                return "vcf"
            return "23andme"
    return "unknown"


def parse_geno_2letter(gt):
    gt = gt.strip().upper()
    alleles = [a for a in gt if a in "ACGT"]
    if not alleles:
        return None
    if len(alleles) == 1:
        alleles *= 2
    return alleles[0], alleles[1]


def run():
    src = config.GENOME_INPUT
    if not src or not src.exists():
        raise SystemExit(f"[step1] нет файла генома. Задай GENOME_INPUT (сейчас: '{src}')")
    fmt = detect_format(src)
    print(f"[step1] вход: {src.name} · формат: {fmt}")

    n_total = n_called = n_nocall = n_amb = 0
    with _open(src) as fin, open(config.NORMALIZED_TYPED, "w", encoding="utf-8") as out:
        out.write("rsid\tchrom\tpos\ta1\ta2\tzygosity\tambiguous\n")
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if fmt == "23andme":
                if len(c) < 4:
                    continue
                rsid, chrom, pos, gt = c[0], c[1].replace("chr", ""), c[2], c[3]
                n_total += 1
                pa = parse_geno_2letter(gt)
                if pa is None:
                    n_nocall += 1
                    continue
                a1, a2 = pa
            else:  # vcf
                if len(c) < 10:
                    continue
                chrom, pos, rsid, ref, alt = c[0].replace("chr", ""), c[1], c[2], c[3], c[4]
                if rsid in (".", ""):
                    rsid = f"{chrom}:{pos}"
                gt = c[9].split(":")[0].replace("|", "/")
                if "." in gt:
                    n_total += 1
                    n_nocall += 1
                    continue
                n_total += 1
                alts = alt.split(",")
                idx = gt.split("/")

                def al(i):
                    if i == "0":
                        return ref
                    return alts[int(i) - 1] if i.isdigit() and int(i) - 1 < len(alts) else None
                a1 = al(idx[0])
                a2 = al(idx[1] if len(idx) > 1 else idx[0])
                if a1 is None or a2 is None or len(a1) != 1 or len(a2) != 1:
                    continue
            zyg = "hom" if a1 == a2 else "het"
            amb = 1 if frozenset({a1, a2}) in config.AMBIGUOUS_PAIRS else 0
            out.write(f"{rsid}\t{chrom}\t{pos}\t{a1}\t{a2}\t{zyg}\t{amb}\n")
            n_called += 1
            n_amb += amb

    print(f"[step1] всего:        {n_total:,}")
    print(f"[step1] валидных:     {n_called:,}")
    print(f"[step1] no-call:      {n_nocall:,}")
    print(f"[step1] палиндромных: {n_amb:,}")
    print(f"[step1] → {config.NORMALIZED_TYPED}")
    return n_called


if __name__ == "__main__":
    run()
