"""Извлечение MHC-региона (chr6) для HIBAG из импутированных dose-файлов, с QC по R².

Правка по ревью (Codex круг 4): раньше mhc_chr6.tsv делался вручную без R²-фильтра и не
регенерировался пайплайном (риск stale + слабая импутация → ложно-уверенные HLA-типы).
Теперь — часть пайплайна, с фильтром R²≥0.5 (для HLA строже) + TYPED всегда.
"""
import gzip
from pathlib import Path
import config

IMP = config.BASE / ".cache" / "imputed"
OUT = config.BASE / ".cache" / "mhc_chr6.tsv"
R2_MIN = 0.5   # для HLA строже, чем общий 0.3


def run():
    fp = IMP / "chr6.dose.vcf.gz"
    if not fp.exists():
        print("[mhc] нет chr6.dose.vcf.gz — пропуск (HIBAG не будет)")
        return False
    n = dropped = 0
    with gzip.open(fp, "rt") as f, open(OUT, "w") as out:
        out.write("pos\tref\talt\tgeno\tr2\n")
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            pos = int(c[1])
            if pos < 28_000_000 or pos > 34_000_000:
                continue
            ref, alt = c[3], c[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            info = c[7]
            typed = "TYPED" in info
            r2 = 0.0
            for fld in info.split(";"):
                if fld.startswith("R2="):
                    try:
                        r2 = float(fld[3:])
                    except ValueError:
                        pass
            if not typed and r2 < R2_MIN:      # QC: слабую импутацию не пускаем в HLA
                dropped += 1
                continue
            s = c[9].split(":")
            gt = s[0].replace("|", "/")
            if gt in (".", "./."):
                continue
            geno = sum(1 for a in gt.split("/") if a == "1")
            out.write(f"{pos}\t{ref}\t{alt}\t{geno}\t{r2 if not typed else 1.0:.3f}\n")
            n += 1
    print(f"[mhc] MHC chr6 SNP (R²≥{R2_MIN} или typed): {n:,} (отброшено по R²: {dropped:,}) → {OUT.name}")
    return True


if __name__ == "__main__":
    run()
