"""ШАГ 13 — Носительство + actionable-находки.

Из значимых ClinVar-находок (P/LP, шаг 5) выделяем:
  (a) гены ACMG SF v3.x — «вторичные находки», по которым есть медицинское действие;
  (b) гетерозиготное носительство P/LP — важно при планировании семьи.

ВАЖНО (честно): SNP-чип покрывает в основном частые варианты, а патогенные ACMG —
редкие, их на чипе обычно нет. Поэтому пустой/короткий результат здесь — нормальный
и ожидаемый исход, а не «всё хорошо проверено по 99 генам». Полнота — только на WGS.
"""
import config

# ACMG SF v3.x — actionable-гены (курированный список ~90 символов)
ACMG_SF = {
    # онко
    "APC", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDKN2A", "EPCAM", "MAX",
    "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH", "NF2", "PALB2", "PMS2", "PTEN",
    "RAD51C", "RAD51D", "RB1", "RET", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD4",
    "STK11", "TMEM127", "TP53", "TSC1", "TSC2", "VHL", "WT1",
    # кардио
    "ACTA2", "ACTC1", "APOB", "BAG3", "CASQ2", "COL3A1", "DES", "DSC2", "DSG2",
    "DSP", "FBN1", "FLNC", "GLA", "KCNH2", "KCNQ1", "LDLR", "LMNA", "MYBPC3",
    "MYH7", "MYH11", "MYL2", "MYL3", "PCSK9", "PKP2", "PRKAG2", "RBM20", "RYR2",
    "SCN5A", "SMAD3", "TGFBR1", "TGFBR2", "TMEM43", "TNNC1", "TNNI3", "TNNT2",
    "TPM1", "TTN", "TTR", "TRDN", "CALM1", "CALM2", "CALM3", "ACVRL1", "ENG",
    # метаболизм/прочее
    "ATP7B", "BTD", "GAA", "OTC", "RYR1", "CACNA1S", "HNF1A", "RPE65", "HFE",
    "FBN2", "MAX",
}


def run():
    rows = []
    with open(config.SIGNIFICANT, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            rec = dict(zip(header, line.rstrip("\n").split("\t")))
            if rec["category"] != "P/LP":
                continue
            rows.append(rec)

    acmg_hits = [r for r in rows if r["gene"] in ACMG_SF]
    carrier_het = [r for r in rows if r["zygosity"] == "het"]

    with open(config.CARRIER, "w", encoding="utf-8") as out:
        out.write("type\tgene\trsid\tclnsig\tstars\tzygosity\tcondition\tmy_a1\tmy_a2\tref\talt\tambiguous\n")
        for r in acmg_hits:
            out.write("ACMG_SF\t" + "\t".join(r.get(c, "") for c in
                      ["gene", "rsid", "clnsig", "stars", "zygosity", "condition",
                       "my_a1", "my_a2", "ref", "alt", "ambiguous"]) + "\n")
        for r in carrier_het:
            out.write("carrier_het\t" + "\t".join(r.get(c, "") for c in
                      ["gene", "rsid", "clnsig", "stars", "zygosity", "condition",
                       "my_a1", "my_a2", "ref", "alt", "ambiguous"]) + "\n")

    print(f"[step13] P/LP всего (несёт аллель): {len(rows)}")
    print(f"[step13]   в actionable-генах ACMG SF: {len(acmg_hits)}")
    print(f"[step13]   гетерозиготное носительство P/LP: {len(carrier_het)}")
    print(f"[step13] → {config.CARRIER}")
    return acmg_hits, carrier_het


if __name__ == "__main__":
    run()
