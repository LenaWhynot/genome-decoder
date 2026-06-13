"""Подготовка чипа к загрузке на Michigan Imputation Server (шаг 2 — импутация).

Берёт Genotek/23andMe VCF (hg19), приводит к формату, который принимают GRCh37-панели
Michigan (HRC r1.1 / 1000G Phase 3): убирает chr-префикс, оставляет аутосомы 1–22 + X,
режет по хромосомам, сортирует, bgzip + tabix. Выход в .cache/michigan-input/ (gitignored).

Сама заливка — вручную пользователем (нужен аккаунт imputationserver.sph.umich.edu + согласие
на отправку генома). Инструкция рядом: MICHIGAN-INSTRUCTIONS.md.
"""
from pathlib import Path
import pysam

BASE = Path(__file__).resolve().parent
SRC = config.GENOME_INPUT   # геном пользователя (raw VCF/txt; для импутации нужен VCF)
OUT = BASE / ".cache" / "michigan-input"
OUT.mkdir(parents=True, exist_ok=True)

# Michigan требует ≥5 образцов → дублируем геном образца в 5 колонок (имутация per-sample
# независима; берём одну колонку из результата). Стандартный обход min-sample гейта.
DUP = 5
SAMPLES = [config.SAMPLE_ID] + [f"{config.SAMPLE_ID}_{i}" for i in range(2, DUP + 1)]
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES) + "\n"
)
KEEP = [str(c) for c in range(1, 23)] + ["X"]   # аутосомы + X (M исключаем)


def run():
    # собрать строки по хромосомам, убрав chr-префикс
    by_chr = {c: [] for c in KEEP}
    n_in = n_kept = 0
    with open(SRC, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0].replace("chr", "")
            if chrom == "M" or chrom == "MT" or chrom not in by_chr:
                continue
            cols[0] = chrom
            by_chr[chrom].append(cols)
            n_kept += 1

    made = []
    for chrom in KEEP:
        rows = by_chr[chrom]
        if not rows:
            continue
        rows.sort(key=lambda c: int(c[1]))          # сортировка по позиции
        plain = OUT / f"chr{chrom}.vcf"
        with open(plain, "w", encoding="utf-8") as out:
            out.write(HEADER)
            for c in rows:
                gt = c[9] if len(c) > 9 else "./."     # генотип образца
                out.write("\t".join(c[:10]) + ("\t" + gt) * (DUP - 1) + "\n")
        gz = str(plain) + ".gz"
        pysam.tabix_compress(str(plain), gz, force=True)   # bgzip
        pysam.tabix_index(gz, preset="vcf", force=True)    # .tbi
        plain.unlink()                                      # оставляем только .gz/.tbi
        made.append((chrom, len(rows)))

    print(f"[michigan] прочитано вариантов: {n_in:,}")
    print(f"[michigan] оставлено (1-22,X, без chr/M): {n_kept:,}")
    print(f"[michigan] файлов по хромосомам: {len(made)}")
    for chrom, n in made:
        print(f"[michigan]   chr{chrom}.vcf.gz : {n:,}")
    print(f"[michigan] → {OUT}")


if __name__ == "__main__":
    run()
