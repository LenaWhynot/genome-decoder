"""ШАГ 7 — Фармакогеномика (PGx).

Аналог PharmCAT: для ключевых генов-ферментов печени определяем по чипу генотип
и фенотип метаболизатора (Poor/Intermediate/Normal/Rapid) + что показать врачу.
Встроенная CPIC-таблица (как embedded-таблица нутригеномики у Григория), plus strand.

Риск-аллели сверены с уже валидированным genome-annotation.md (LPA-урок: не «на глаз»).
Это не назначение — карточка для врача.
"""
import config

# rsid -> (ген, аллель варианта на plus strand, метка стар-аллеля/эффект, препараты)
PGX = {
    # CYP2C19 — клопидогрел, ИПП, СИОЗС
    "rs4244285":   ("CYP2C19", "A", "*2 LoF", "клопидогрел, ИПП, СИОЗС"),
    "rs4986893":   ("CYP2C19", "A", "*3 LoF", "клопидогрел, ИПП, СИОЗС"),
    "rs12248560":  ("CYP2C19", "T", "*17 gain", "клопидогрел, ИПП, СИОЗС"),
    # CYP2C9 — НПВС, варфарин, лозартан, фенитоин
    "rs1799853":   ("CYP2C9", "T", "*2 reduced", "НПВС, варфарин, лозартан"),
    "rs1057910":   ("CYP2C9", "C", "*3 reduced", "НПВС, варфарин, лозартан, фенитоин"),
    # VKORC1 — чувствительность к варфарину
    "rs9923231":   ("VKORC1", "A", "пром. сниж. дозы", "варфарин"),
    # SLCO1B1 — миопатия на статинах
    "rs4149056":   ("SLCO1B1", "C", "*5 decreased", "симвастатин и др. статины"),
    # TPMT — тиопурины
    "rs1142345":   ("TPMT", "C", "*3C", "азатиоприн, 6-меркаптопурин"),
    "rs1800460":   ("TPMT", "A", "*3B", "азатиоприн, 6-меркаптопурин"),
    # NUDT15 — тиопурины (важно для не-европейцев, но проверяем)
    "rs116855232": ("NUDT15", "T", "*3 LoF", "азатиоприн, 6-меркаптопурин"),
    # DPYD — токсичность 5-FU/капецитабина
    "rs3918290":   ("DPYD", "A", "*2A LoF", "5-фторурацил, капецитабин"),
    "rs55886062":  ("DPYD", "C", "*13 LoF", "5-фторурацил, капецитабин"),
    "rs67376798":  ("DPYD", "A", "c.2846 reduced", "5-фторурацил, капецитабин"),
    # CYP3A5 — такролимус
    "rs776746":    ("CYP3A5", "A", "*1 expressor", "такролимус"),
    # CYP1A2 — кофеин/тео­филлин (медленный индуктор)
    "rs762551":    ("CYP1A2", "A", "*1F быстрый", "кофеин, теофиллин"),
    # UGT1A1 — иринотекан/атазанавир + синдром Жильбера
    "rs887829":    ("UGT1A1", "T", "*80 сниж.", "иринотекан, атазанавир"),
}

# простые фенотип-каллы по числу LoF-копий
def phenotype(gene, calls):
    """calls: {rsid: alt_copies}. Возвращает строку фенотипа или ''."""
    if gene == "CYP2C9":
        red = calls.get("rs1799853", 0) + calls.get("rs1057910", 0)
        return {0: "Normal metabolizer", 1: "Intermediate metabolizer",
                2: "Poor metabolizer"}.get(min(red, 2), "")
    if gene == "CYP2C19":
        lof = calls.get("rs4244285", 0) + calls.get("rs4986893", 0)
        gain = calls.get("rs12248560", 0)
        if lof >= 2:
            return "Poor metabolizer"
        if lof == 1 and gain == 0:
            return "Intermediate metabolizer"
        if lof == 0 and gain >= 1:
            return "Rapid/Ultrarapid metabolizer"
        if lof == 0 and gain == 0:
            return "Normal metabolizer"
        return "Intermediate (смешанный *17)"
    if gene == "SLCO1B1":
        c = calls.get("rs4149056", 0)
        return {0: "Normal function", 1: "Decreased (риск миопатии↑)",
                2: "Poor function (риск миопатии↑↑)"}.get(c, "")
    if gene == "CYP1A2":
        # rs762551: A=*1F (быстрый индуцируемый), C/C=медленный
        a = calls.get("rs762551", 0)
        return {2: "*1F/*1F быстрый индуцируемый", 1: "промежуточный",
                0: "медленный (C/C)"}.get(a, "")
    return ""


def load_norm():
    d = {}
    with open(config.NORMALIZED_TYPED, encoding="utf-8") as f:
        next(f)
        for line in f:
            rsid, chrom, pos, a1, a2, zyg, amb = line.rstrip("\n").split("\t")
            d[rsid] = (a1, a2)
    return d


def run():
    norm = load_norm()
    by_gene = {}
    rows = []
    for rsid, (gene, var, label, drugs) in PGX.items():
        gt = norm.get(rsid)
        if gt is None:
            rows.append((gene, rsid, "—", "нет на чипе V5", var, label, drugs))
            continue
        a1, a2 = gt
        copies = [a1, a2].count(var)
        by_gene.setdefault(gene, {})[rsid] = copies
        carry = {0: "нет", 1: "гетерозигота", 2: "гомозигота"}[copies]
        rows.append((gene, rsid, f"{a1}{a2}", carry, var, label, drugs))

    with open(config.PGX, "w", encoding="utf-8") as out:
        out.write("gene\trsid\tgenotype\tcarries_variant\tvariant_allele\tstar_effect\tdrugs\tgene_phenotype\n")
        for gene, rsid, gtxt, carry, var, label, drugs in rows:
            ph = phenotype(gene, by_gene.get(gene, {}))
            out.write(f"{gene}\t{rsid}\t{gtxt}\t{carry}\t{var}\t{label}\t{drugs}\t{ph}\n")

    print(f"[step7] PGx-маркеров проверено: {len(PGX)}")
    print("[step7] фенотипы:")
    for gene in ["CYP2C9", "CYP2C19", "SLCO1B1", "CYP1A2"]:
        ph = phenotype(gene, by_gene.get(gene, {}))
        if ph:
            print(f"[step7]   {gene}: {ph}")
    print(f"[step7] → {config.PGX}")
    return rows


if __name__ == "__main__":
    run()
