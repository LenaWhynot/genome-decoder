"""ШАГ 12 — HLA-типирование (таргетная тег-SNP панель).

Не полный HIBAG, а прагматичная панель: клинически важные HLA-аллели через хорошо
валидированные тег-SNP. Закрывает actionable-часть — БЕЗОПАСНОСТЬ ЛЕКАРСТВ (тяжёлые
кожные реакции) + целиакия. Читает типизированные данные по rsID (01_normalized_typed.tsv).

⚠️ Тег-SNP = вероятностная привязка, не прямое HLA-типирование. Положительный тег →
подтверждать клиническим HLA-тестом ПЕРЕД назначением препарата. Часть тегов
популяционно-специфична (B*15:02 — надёжен у восточноазиатов).
"""
import config

TYPED = config.OUT / "01_normalized_typed.tsv"

# rsid -> (HLA-аллель, тег-аллель, препарат/состояние, риск при носительстве)
HLA = {
    "rs2395029":  ("HLA-B*57:01", "G", "абакавир (ВИЧ), флуклоксациллин",
                   "тяжёлая гиперчувствительность — препарат НЕ назначать"),
    "rs9263726":  ("HLA-B*58:01", "A", "аллопуринол (подагра)",
                   "риск SJS/TEN (тяжёлые кожные реакции)"),
    # HLA-A*31:01 (карбамазепин) НЕ включён: надёжного однонуклеотидного тега нет
    # (rs1061235 = частый аллель → переоценивает). Нужно прямое HLA-типирование.
    "rs3909184":  ("HLA-B*15:02", "A", "карбамазепин/фенитоин (вост.-азиаты)",
                   "риск SJS/TEN — тег надёжен в азиатских популяциях"),
    "rs2187668":  ("HLA-DQ2.5", "T", "глютен (целиакия)",
                   "предрасположенность к целиакии (DQA1*05:01-DQB1*02:01)"),
    "rs7454108":  ("HLA-DQ8", "C", "глютен (целиакия)",
                   "предрасположенность к целиакии (DQB1*03:02)"),
}


def load_typed():
    d = {}
    try:
        with open(TYPED, encoding="utf-8") as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                d[p[0]] = (p[3], p[4])      # rsid -> (a1,a2)
    except FileNotFoundError:
        print(f"[step12] нет {TYPED} — нужен типизированный нормализованный (бэкап)")
    return d


def run():
    typed = load_typed()
    rows = []
    for rsid, (allele, tag, drug, risk) in HLA.items():
        gt = typed.get(rsid)
        if gt is None:
            rows.append((allele, rsid, "—", "нет на чипе", drug, risk))
            continue
        a1, a2 = gt
        copies = [a1, a2].count(tag)
        carry = {0: "тег отрицателен", 1: "тег+ (гетерозигота)", 2: "тег+ (гомозигота)"}[copies]
        rows.append((allele, rsid, f"{a1}{a2}", carry, drug, risk))

    with open(config.OUT / "12_hla.tsv", "w", encoding="utf-8") as out:
        out.write("hla_allele\trsid\tgenotype\ttag_status\tdrug\trisk\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    pos = [r for r in rows if "тег+" in r[3]]
    print(f"[step12] HLA тег-SNP проверено: {len(HLA)}")
    for r in rows:
        mark = "⚠" if "тег+" in r[3] else " "
        print(f"[step12]  {mark} {r[0]}: {r[3]} ({r[2]}) — {r[4]}")
    print(f"[step12] положительных тегов: {len(pos)} → out/12_hla.tsv")
    return rows


if __name__ == "__main__":
    run()
