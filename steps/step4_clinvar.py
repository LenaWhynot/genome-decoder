"""ШАГ 4 — ClinVar.

Сверяем ВСЕ типизированные SNP (не 37 вручную) с мировой базой клинических оценок.
Скачиваем тиражируемую ClinVar GRCh37 VCF (build 37 = совпадает с 23andMe), стримим,
для каждого совпадения по chrom:pos определяем — несёт ли образец клинически оценённый
аллель (ALT), и в какой зиготности. Сохраняем оценку + уровень доверия (звёзды).
"""
import gzip
import urllib.request
import config

# CLNREVSTAT → «золотые звёзды» доверия ClinVar
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
}


def download_clinvar():
    if config.CLINVAR_VCF_GZ.exists() and config.CLINVAR_VCF_GZ.stat().st_size > 1_000_000:
        mb = config.CLINVAR_VCF_GZ.stat().st_size / 1e6
        print(f"[step4] ClinVar уже скачан ({mb:.0f} MB), пропускаю загрузку")
        return
    print(f"[step4] скачиваю ClinVar GRCh37 → {config.CLINVAR_VCF_GZ} ...")
    urllib.request.urlretrieve(config.CLINVAR_URL, config.CLINVAR_VCF_GZ)
    mb = config.CLINVAR_VCF_GZ.stat().st_size / 1e6
    print(f"[step4] готово: {mb:.0f} MB")


def load_my_variants():
    """{(chrom,pos): [ (rsid,a1,a2,zyg,amb,r2,typed), ... ]} — СПИСОК записей на позицию.
    Список, а не одна запись — чтобы multiallelic-сайты (несколько ALT на одной позиции)
    не перетирали друг друга. r2/typed опциональны (есть в импутированном)."""
    idx = {}
    with open(config.NORMALIZED, encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rsid, chrom, pos, a1, a2, zyg, amb = p[:7]
            r2 = p[7] if len(p) > 7 else ""
            typed = p[8] if len(p) > 8 else ""
            idx.setdefault((chrom, pos), []).append((rsid, a1, a2, zyg, amb, r2, typed))
    return idx


def parse_info(info: str):
    d = {}
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
    return d


def stars_for(revstat: str) -> int:
    return STARS.get(revstat, 0)


def max_af(inf: dict):
    """Максимальная популяционная частота из зашитых в ClinVar полей (ESP/ExAC/1000G).
    Это частотный фильтр здравого смысла (шаг 5) без скачивания полной gnomAD."""
    best = ""
    for k in ("AF_EXAC", "AF_TGP", "AF_ESP"):
        v = inf.get(k)
        if v:
            try:
                fv = float(v)
                if best == "" or fv > best:
                    best = fv
            except ValueError:
                pass
    return "" if best == "" else f"{best:.5f}"


def run():
    download_clinvar()
    mine = load_my_variants()
    print(f"[step4] моих типизированных позиций: {len(mine):,}")

    n_records = n_pos_match = n_carry = 0
    with gzip.open(config.CLINVAR_VCF_GZ, "rt", encoding="utf-8", errors="replace") as f, \
         open(config.CLINVAR_HITS, "w", encoding="utf-8") as out:
        out.write("rsid\tchrom\tpos\tref\talt\tmy_a1\tmy_a2\talt_copies\t"
                  "zygosity\tambiguous\tclnsig\tstars\tmax_af\tgene\tcondition\tclnvc\tr2\ttyped\n")
        for line in f:
            if line.startswith("#"):
                continue
            n_records += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, _q, _flt, info = cols[:8]
            key = (chrom, pos)
            hit_list = mine.get(key)
            if not hit_list:
                continue
            n_pos_match += 1
            # КРИТИЧНО: чип/импутация надёжно сравниваются только по SNV. Инделы/MNV (ref/alt
            # длиннее 1 буквы) — НЕЛЬЗЯ (давало ложные «гомозиготы» ref=TA alt=T vs gt=TT).
            if not (len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT"):
                continue
            # multiallelic: ищем СРЕДИ её записей на этой позиции ту, что несёт именно этот ALT
            matched = None
            for rec in hit_list:
                if alt in (rec[1], rec[2]):     # a1 или a2 == ClinVar ALT
                    matched = rec
                    break
            if matched is None:
                continue  # она не несёт оценённый (ALT) аллель этой ClinVar-записи
            rsid, a1, a2, zyg, amb, r2, typed = matched
            alt_copies = [a1, a2].count(alt)
            n_carry += 1
            inf = parse_info(info)
            clnsig = inf.get("CLNSIG", inf.get("CLNSIGINCL", "")).replace("_", " ")
            revstat = inf.get("CLNREVSTAT", "")
            stars = stars_for(revstat)
            gene = inf.get("GENEINFO", "").split(":")[0].split("|")[0]
            cond = inf.get("CLNDN", "").replace("_", " ")[:120]
            clnvc = inf.get("CLNVC", "")
            af = max_af(inf)
            zyg_alt = "hom" if alt_copies == 2 else "het"
            out.write(f"{rsid}\t{chrom}\t{pos}\t{ref}\t{alt}\t{a1}\t{a2}\t{alt_copies}\t"
                      f"{zyg_alt}\t{amb}\t{clnsig}\t{stars}\t{af}\t{gene}\t{cond}\t{clnvc}\t{r2}\t{typed}\n")

    print(f"[step4] записей ClinVar просмотрено: {n_records:,}")
    print(f"[step4] совпало по позиции:          {n_pos_match:,}")
    print(f"[step4] из них образец несёт ALT-аллель: {n_carry:,}")
    print(f"[step4] → {config.CLINVAR_HITS}")
    return n_carry


if __name__ == "__main__":
    run()
