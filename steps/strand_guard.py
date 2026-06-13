"""Strand guard для рукописных effect-аллелей (шаги 9/10).

Проблема. В embedded-таблицах трейтов/нутригеномики effect-букву легко
записать на ЧУЖОЙ цепи (минус / mRNA-sense), тогда как генотипы из чипа/VCF —
на плюс-цепи GRCh37. При расхождении `count(eff)` молча даёт 0 копий, и текст
подбирается ЗЕРКАЛЬНО — у мажорного гомозиготы это маскируется под «верный»
ответ. Так ловятся флипы вида CLOCK eff="C" или GC eff="C" (нужно "G"), а у
HERC2/OCA2 (цвет глаз) перевёрнутый словарь даёт «карие» там, где должно быть
«голубые».

Решение. Для КАЖДОГО рукописного effect-аллеля шага 9/10 проверяем, что он
принадлежит множеству ПЛЮС-ЦЕПНЫХ аллелей позиции (Ensembl forward GRCh37).
Если буквы нет, а её комплемент есть — это strand-флип. По умолчанию падаем
жёстко; с STRAND_GUARD_WARN=1 — печатаем видимое предупреждение и продолжаем.

При добавлении нового rsID в шаги 9/10 впишите его плюс-цепные аллели в
PLUS_STRAND_ALLELES (берётся из Ensembl REST allele_string, forward strand) —
страж потребует этого явно, а не «на глаз».
"""
import os
import sys

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# rsID -> множество ПЛЮС-ЦЕПНЫХ аллелей позиции (Ensembl forward GRCh37).
PLUS_STRAND_ALLELES = {
    # ── шаг 9 (трейты) ──
    "rs1815739":  {"C", "T"},        # ACTN3
    "rs4680":     {"G", "A"},        # COMT
    "rs6265":     {"C", "T"},        # BDNF
    "rs1801260":  {"A", "G"},        # CLOCK  (буквы "C" тут нет — ловит флип)
    "rs12913832": {"A", "G"},        # HERC2/OCA2
    "rs1426654":  {"A", "G"},        # SLC24A5
    "rs17822931": {"C", "T"},        # ABCC11
    "rs53576":    {"A", "G", "T"},   # OXTR
    # ── шаг 10 (нутригеномика) ──
    "rs4988235":  {"G", "A"},        # LCT/MCM6
    "rs762551":   {"C", "A"},        # CYP1A2
    "rs1229984":  {"C", "T"},        # ADH1B
    "rs671":      {"G", "A"},        # ALDH2
    "rs1801133":  {"G", "A"},        # MTHFR C677T
    "rs174537":   {"G", "C", "T"},   # FADS1
    "rs9939609":  {"T", "A"},        # FTO
    "rs7903146":  {"C", "T"},        # TCF7L2
    "rs12934922": {"A", "T"},        # BCO1
    "rs1800562":  {"G", "A"},        # HFE
    "rs2282679":  {"T", "G"},        # GC  (буквы "C" тут нет — ловит флип)
    "rs713598":   {"C", "G", "T"},   # TAS2R38
    "rs2187668":  {"C", "T"},        # HLA-DQ2.5 tag
}


def check_effect_alleles(table, norm, step_name):
    """Проверить, что каждый рукописный eff лежит на плюс-цепи позиции.

    table : dict rsid -> (..., eff, ...) — eff на индексе 2 (как в TRAITS/NUTRI).
    norm  : dict rsid -> (a1, a2) из нормализованного файла (наблюдённые плюс-цепные).
    Возвращает список строк-флагов (для печати), пустой если чисто.
    Падает (SystemExit), если найден strand-флип и STRAND_GUARD_WARN != "1".
    """
    warn_only = os.environ.get("STRAND_GUARD_WARN") == "1"
    flags = []
    for rsid, spec in table.items():
        eff = spec[2]
        allowed = PLUS_STRAND_ALLELES.get(rsid)

        if allowed is None:
            msg = (f"[{step_name}] СТРАЖ: для {rsid} не задано множество плюс-цепных "
                   f"аллелей в strand_guard.PLUS_STRAND_ALLELES — впиши перед использованием")
            flags.append(msg)
            if not warn_only:
                sys.exit("❌ " + msg)
            continue

        if eff in allowed:
            continue  # effect-буква на плюс-цепи — ок

        comp = COMPLEMENT.get(eff)
        if comp in allowed:
            msg = (f"[{step_name}] СТРАЖ: effect-аллель {rsid} eff='{eff}' НЕ на плюс-цепи "
                   f"(там {sorted(allowed)}); похоже на чужую нить — нужен '{comp}'. STRAND-ФЛИП.")
        else:
            msg = (f"[{step_name}] СТРАЖ: effect-аллель {rsid} eff='{eff}' не из плюс-цепных "
                   f"{sorted(allowed)} и комплемент не подходит — проверь rsID/координату.")
        flags.append(msg)
        if not warn_only:
            sys.exit("❌ " + msg)

        gt = norm.get(rsid)
        if gt and eff not in gt and comp in gt:
            flags.append(f"[{step_name}] СТРАЖ: подтверждено генотипом — {rsid} call={gt}, "
                         f"видим '{comp}', не '{eff}'.")

    return flags
