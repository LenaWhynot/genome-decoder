"""ШАГ 11 — Происхождение (прагматично, без 1000G-download).

Не полный PCA/ADMIXTURE, а честный сигнал: панель сильных ancestry-informative маркеров
(AIM) с известным «европейским» аллелем + подсчёт mtDNA-маркеров. Y-гаплогруппы нет
(женщина). Точная mtDNA-гаплогруппа и % примеси требуют Haplogrep/1000G — здесь не
претендуем (риск ошибки), даём направленный вывод.
"""
import config

TYPED = config.OUT / "01_normalized_typed.tsv"
SRC_VCF = config.GENOME_INPUT

# rsid -> (ген, «европейский»/информативный аллель, что значит)
AIM = {
    "rs1426654":  ("SLC24A5", "A", "светлая кожа — A≈100% у европейцев, ≈0% у африканцев"),
    "rs16891982": ("SLC45A2", "G", "светлая кожа — G≈98% у европейцев"),
    "rs3827760":  ("EDAR", "A", "A = НЕ восточноазиатский (G370A — азиатский маркер)"),
    "rs2814778":  ("DARC/Duffy", "T", "T = НЕ африканский (C = Duffy-null, африканский)"),
    "rs1129038":  ("HERC2", "T", "T = голубоглазый/европейский вариант"),
    # rs12913832 убран по ревью (Codex круг 4): это маркер цвета глаз, не ancestry-AIM,
    # + strand-неоднозначность A/G. Пигментация уже в трейтах (step9).
}


def load_typed():
    d = {}
    try:
        with open(TYPED, encoding="utf-8") as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                d[p[0]] = (p[3], p[4])
    except FileNotFoundError:
        pass
    return d


def count_mtdna():
    n = 0
    try:
        with open(SRC_VCF, encoding="utf-8") as f:
            for line in f:
                if line.startswith("chrM\t") or line.startswith("MT\t"):
                    n += 1
    except FileNotFoundError:
        pass
    return n


def run():
    typed = load_typed()
    rows = []
    eur_hits = 0
    for rsid, (gene, eur, note) in AIM.items():
        gt = typed.get(rsid)
        if gt is None:
            rows.append((gene, rsid, "—", "нет на чипе", note))
            continue
        a1, a2 = gt
        copies = [a1, a2].count(eur)
        status = {2: "гомозигота (евр.)", 1: "гетерозигота", 0: "не несёт евр. аллель"}[copies]
        if copies == 2:
            eur_hits += 1
        rows.append((gene, rsid, f"{a1}{a2}", status, note))

    mt = count_mtdna()
    with open(config.OUT / "11_ancestry.tsv", "w", encoding="utf-8") as out:
        out.write("gene\trsid\tgenotype\tstatus\tnote\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
        out.write(f"mtDNA\t—\t{mt} маркеров\tпрофиль есть, гаплогруппа = Haplogrep\t—\n")

    typed_aims = [r for r in rows if r[3] != "нет на чипе"]
    print(f"[step11] AIM проверено: {len(typed_aims)}/{len(AIM)}")
    for r in rows:
        print(f"[step11]   {r[0]} {r[1]}: {r[3]} ({r[2]})")
    print(f"[step11] гомозигот по европ. аллелю: {eur_hits}/{len(typed_aims)} → сильный европейский сигнал")
    print(f"[step11] mtDNA: {mt} маркеров (Y нет — женщина); точная гаплогруппа = Haplogrep")
    print(f"[step11] → out/11_ancestry.tsv")
    return rows


if __name__ == "__main__":
    run()
