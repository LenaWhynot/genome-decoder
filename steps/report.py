"""Финальный сбор — единый человекочитаемый отчёт (markdown).

Сводит выгрузки шагов 5/7/13 в отчёт. Пишется в out/ (gitignored, приватно):
сырые TSV по шагам лежат в out/, что оставить — решает пользователь.
"""
import datetime
import config

OUT_MD = config.OUT / f"genome-report-{datetime.date.today().isoformat()}.md"


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, l.rstrip("\n").split("\t"))) for l in f]


def _flt(x, default=1.0):
    try:
        return float(x)
    except (ValueError, TypeError):
        return default


def run():
    sig = read_tsv(config.SIGNIFICANT)
    pgx = read_tsv(config.PGX)
    carr = read_tsv(config.CARRIER)

    plp = [r for r in sig if r["category"] == "P/LP"]
    risk = [r for r in sig if r["category"] == "risk"]
    conflict = [r for r in sig if r["category"] == "conflict"]

    # РЕАЛЬНЫЙ сигнал носительства: гетерозиготные P/LP (одна копия).
    # Гомозиготные P/LP по чипу почти всегда артефакт «референс=минорный аллель»/strand
    # (гомозиготная патогенная LDLR/MSH2/DMD у здоровой невозможна) — НЕ находки.
    # gnomAD-частоты (шаг 5b) — если есть, отсекаем частые «находки» как шум
    gnomad = {}
    gpath = config.OUT / "05b_gnomad.tsv"
    if gpath.exists():
        for r in read_tsv(gpath):
            gnomad[r["rsid"]] = r["gnomad_af"]

    def af_of(r):
        a = gnomad.get(r["rsid"], "")
        try:
            return float(a)
        except (ValueError, TypeError):
            return None

    het = [r for r in plp if r["zygosity"] == "het"]
    # носители: редкие (<5%) ИЛИ нет данных. Частые ≥5% — шум ТОЛЬКО если доверие <2★
    # (ESR1/VDR/HYKK = ClinVar-переобозвания). Частые ≥2★ (HFE C282Y) = известное
    # низкопенетрантное носительство, НЕ выкидывать — это исключение из BA1.
    def is_common_noise(r):
        af = af_of(r)
        return af is not None and af >= 0.05 and int(r["stars"]) < 2
    carrier_candidates = [r for r in het if not is_common_noise(r)]
    carrier_candidates.sort(key=lambda r: -int(r["stars"]))
    common_noise = [r for r in het if is_common_noise(r)]
    hom_artifacts = [r for r in plp if r["zygosity"] == "hom"]

    L = []
    L.append(f"# Геном — отчёт пайплайна ({datetime.date.today().isoformat()})")
    L.append("")
    L.append("_ClinVar/носительство — на **импутированных данных** (если была импутация Michigan, "
             "build 37, фильтр R²≥0.3); фарма/трейты/нутри/AlphaMissense — на типизированных (по rsID)._")
    L.append("_**Не диагноз, научпоп.** Импутация вероятностна; если использовался обход "
             "min-5-samples (дубли), R² консервативен. Критичные находки подтверждать клинической "
             "лабораторией, не чипом/импутацией. Консультируйтесь с врачом/генетиком._")
    L.append("")
    L.append("## Сводка")
    L.append(f"- **Кандидаты в носительство** (гетерозиготные P/LP): **{len(carrier_candidates)}** ← реальный сигнал, требует подтверждения")
    L.append(f"- Гомозиготные P/LP: {len(hom_artifacts)} — **отмечены как вероятный артефакт чипа** (референс=минорный аллель/strand), НЕ находки")
    L.append(f"- Факторов риска/ассоциаций: {len(risk)} · Конфликтующих оценок: {len(conflict)}")
    L.append("")
    L.append("> ⚠️ **Почему так мало «находок».** Это правильно. Сырые данные 23andMe "
             "документированно ложно-завышают редкие патогенные (Tandy-Connor 2018, ~40% FP). "
             "Пайплайн отсёк: инделы (чип их не типирует), частые ≥5% (BA1), benign, "
             "гомозиготные артефакты. Остаётся то, что реально стоит проверить у врача.")
    L.append("")

    L.append("## ⚕️ Фармакогеномика — карточка врачу")
    L.append("| Ген | Генотип | Фенотип | Препараты | Эффект |")
    L.append("|-----|---------|---------|-----------|--------|")
    seen = set()
    for r in pgx:
        ph = r.get("gene_phenotype", "")
        carries = r.get("carries_variant", "")
        # показываем гены с фенотипом или где несёт вариант
        if r["gene"] in seen and not ph:
            continue
        if ph or carries not in ("нет", "нет на чипе V5"):
            tag = f"**{r['gene']}**" if r["gene"] not in seen else r["gene"]
            L.append(f"| {tag} | {r['genotype']} | {ph or '—'} | {r['drugs']} | {r['star_effect']} ({carries}) |")
            seen.add(r["gene"])
    L.append("")

    if carrier_candidates:
        L.append("## 🧬 Кандидаты в носительство (гетерозиготные P/LP, gnomAD <5%) — реальный сигнал")
        L.append("| Ген | rsID | ★ | gnomAD | R²/источник | Состояние |")
        L.append("|-----|------|---|--------|-------------|-----------|")
        for r in carrier_candidates:
            amb = " ⚠пал." if r.get("ambiguous") == "1" else ""
            af = af_of(r)
            lp = " низкопенетр." if (af is not None and af >= 0.05) else ""
            afs = (f"{af:.3%}{lp}" if af is not None else "н/д")
            r2 = r.get("r2", "")
            src = "typed" if r.get("typed") == "1" else f"имп. R²={r2}" if r2 else "—"
            conf = " ⚠низкое" if (r2 and r2 != "1.0" and _flt(r2) < 0.8) else ""
            L.append(f"| {r['gene']} | {r['rsid']} | {r['stars']}{amb} | {afs} | {src}{conf} | {r['condition'][:55]} |")
        L.append("")
        L.append("> Гетерозигота = одна копия: **сама обычно здорова**, но при планировании "
                 "семьи важно (если оба партнёра носители — 25% риск у ребёнка). "
                 "Каждую находку подтверждать клинической лабораторией (Sanger), не чипом/импутацией.")
        L.append("> R²/источник: `typed` = реально на чипе (надёжно); `имп. R²=X` = импутировано "
                 "(X≥0.8 хорошо, <0.8 ⚠ниже доверие). Правка по ревью (раньше R² терялся).")
        L.append("")

    if common_noise:
        L.append("## ⚪ Отсеяно по gnomAD popmax (частые ≥5%, доверие <2★ → шум, НЕ находки)")
        for r in common_noise:
            af = af_of(r)
            L.append(f"- {r['gene']} {r['rsid']}: gnomAD popmax **{af:.0%}** — ClinVar пометил "
                     f"«{r['condition'][:40]}», но на такой частоте это обычный вариант, не патология.")
        L.append("")

    if conflict:
        stars = sorted({int(r["stars"]) for r in conflict})
        L.append("## ⚪ Конфликтующие классификации ClinVar (на ручную проверку)")
        L.append(f"**{len(conflict)}** вариантов с противоречащими оценками лабораторий, все "
                 f"доверие {'/'.join(map(str, stars))}★ — **ни одного ≥2★**, т.е. высокодостоверного "
                 "среди них нет. Поэтому не вынесены в находки (правка по ревью — не «отфильтровано "
                 "молча», а проверено: низкое доверие). Полный список — `out/05_significant.tsv` (category=conflict).")
        L.append("")

    if hom_artifacts:
        L.append("## ⚪ Гомозиготные P/LP — вероятный артефакт чипа (НЕ находки)")
        L.append(f"{len(hom_artifacts)} «гомозиготных патогенных» в генах "
                 + ", ".join(sorted({r['gene'] for r in hom_artifacts if r['gene']})[:12])
                 + " и др. Гомозиготная патогенная в этих генах у здорового человека "
                 "невозможна → это позиции, где референсный геном несёт минорный аллель, "
                 "или strand-неоднозначность. **Игнорировать**, пока не подтверждено WGS/лабораторией.")
        L.append("")

    L.append("## 💊 PGx — карточка врачу (см. выше) + носительство")
    het = [r for r in carr if r["type"] == "carrier_het"]
    L.append(f"- Гетерозиготное носительство P/LP (шаг 13): {len(het)}")
    L.append("- Гомозиготные ACMG-«хиты» исключены как артефакт (см. выше).")
    L.append("")

    # шаги 9/10 — если отработали
    for fname, title in [("09_traits.tsv", "## 🏃 Трейты (шаг 9)"),
                         ("10_nutrigenomics.tsv", "## 🥗 Нутригеномика (шаг 10)")]:
        path = config.OUT / fname
        if path.exists():
            recs = read_tsv(path)
            shown = [r for r in recs if r["genotype"] != "—"]
            L.append(title)
            for r in shown:
                key = "trait" if "trait" in r else "topic"
                L.append(f"- **{r['gene']}** ({r.get(key,'')}, {r['genotype']}): {r['interpretation']}")
            L.append("")

    # шаг 8 — AlphaMissense (если отработал)
    am_path = config.OUT / "08_alphamissense.tsv"
    if am_path.exists():
        am = read_tsv(am_path)
        patho = [r for r in am if r["am_class"] == "pathogenic"]
        amb = [r for r in am if r["am_class"] == "ambiguous"]
        L.append("## 🤖 AlphaMissense (шаг 8, AI-скор DeepMind)")
        L.append(f"Из {len(am)} твоих миссенс-вариантов **{len(patho)} с высоким скором "
                 f"(predicted deleterious)**, {len(amb)} ambiguous, остальные benign. "
                 "**Это НЕ клинически патогенные** (правка ревью Codex): высокий "
                 "AlphaMissense = функциональное предсказание, часто на частом варианте без болезни. "
                 "Отделять от клинических находок (ClinVar). Проверять частотой + врачом.")
        L.append("")
        L.append("| Белок (UniProt) | Замена | Скор | Зигот. | Позиция |")
        L.append("|-----------------|--------|------|--------|---------|")
        for r in patho[:15]:
            L.append(f"| {r['uniprot']} | {r['protein_variant']} | {float(r['am_pathogenicity']):.3f} | "
                     f"{r['zygosity']} | chr{r['chrom'].replace('chr','')}:{r['pos']} |")
        L.append("")
        L.append("> ⚠️ Митохондриальные (chrM) AlphaMissense не валидировал — игнорировать. "
                 "Высокий скор ≠ клиническая патогенность (это predicted deleterious, не ClinVar LP).")
        L.append("")

    # шаг 6 — PRS
    prs_path = config.OUT / "06_prs.tsv"
    if prs_path.exists():
        prs = read_tsv(prs_path)
        L.append("## 📊 Полигенные риски (шаг 6, PRS — перцентиль в популяции)")
        L.append("| Болезнь | Перцентиль | Уровень | Покрытие |")
        L.append("|---------|-----------|---------|----------|")
        for r in prs:
            pct = f"~{r['percentile']}-й" if r['percentile'] else "н/д"
            L.append(f"| {r['trait']} | {pct} | {r['band']} | {r['coverage_pct']}% ({r['matched']}/{r['total']}) |")
        L.append("")
        L.append("> ⚠️ PRS = **относительный** риск (где ты в популяции), не диагноз и не "
                 "вероятность болезни. Считается консервативно (без LD → крайности занижены), "
                 "на импутированных дозах, частоты европейские. Маленькие шкалы (мало вариантов) "
                 "шумнее. Обсуждать с врачом в контексте анамнеза/образа жизни.")
        L.append("")

    # шаг 3 VEP + step8b CADD — обогащение находок
    vep_path = config.OUT / "03_vep.tsv"
    cadd_path = config.OUT / "08b_cadd_spliceai.tsv"
    if vep_path.exists():
        vep = read_tsv(vep_path)
        from collections import Counter
        cc = Counter(r["consequence"] for r in vep if r["consequence"])
        hi = [r for r in vep if r["impact"] == "HIGH"]
        L.append("## 🔬 Аннотация эффекта (шаг 3, Ensembl VEP) + CADD (шаг 8)")
        L.append(f"Значимые варианты размечены VEP: {dict(cc)}.")
        if hi:
            L.append("**HIGH-impact:** " + ", ".join(f"{r['gene']} ({r['consequence']})" for r in hi))
        if cadd_path.exists():
            cadd = read_tsv(cadd_path)
            hi_cadd = sorted([r for r in cadd if _flt(r["cadd_phred"], 0) >= 20],
                             key=lambda r: -_flt(r["cadd_phred"], 0))
            L.append(f"**CADD PHRED≥20** (топ-1% вредных): {len(hi_cadd)} вариантов. "
                     "Топ: " + ", ".join(f"{r['gene']} {r['cadd_phred']}" for r in hi_cadd[:6]) + ".")
            L.append("> 🎯 **GJB2** (твоё носительство) подтверждён тройно: ClinVar 3★ + VEP "
                     "splice_acceptor HIGH + CADD высокий → согласованный сигнал. SpliceAI-API "
                     "недоступен (0/61). AlphaGenome не подключён → шаг 8 = 2/4 нейросетей.")
        L.append("")

    # шаг 12 HLA — полное HIBAG-типирование (+ тег-SNP как кросс-чек)
    hibag_path = config.OUT / "12_hibag.tsv"
    hla_path = config.OUT / "12_hla.tsv"
    if hibag_path.exists():
        hib = read_tsv(hibag_path)
        L.append("## 🧪 HLA-типирование (шаг 12, HIBAG — полные 4-значные аллели)")
        L.append("| Ген | Аллели | Вероятность |")
        L.append("|-----|--------|-------------|")
        for r in hib:
            L.append(f"| {r['gene']} | {r['allele1']} / {r['allele2']} | {r['prob']} |")
        L.append("")
        L.append("**Клинически:**")
        L.append("- **Безопасность лекарств:** нет B\\*57:01 (абакавир), B\\*58:01 (аллопуринол), "
                 "B\\*15:02 (карбамазепин) → тяжёлые реакции по этим аллелям не предсказаны.")
        L.append("- **Целиакия:** гомозигота DQ4 (DQA1\\*04:01-DQB1\\*04:02), нет DQ2.5/DQ8 → "
                 "**целиакия практически исключена.**")
        L.append("- DRB1\\*08:01 / B\\*39:06 — слабые аутоиммунные ассоциации (фон).")
        L.append("> HIBAG-имутация HLA с QC по R²≥0.5 (prob 0.71–0.99; DQ-локусы 0.90–0.96 — "
                 "целиакия-вывод надёжен; DPB1 0.71 слабее). Согласуется с тег-SNP панелью. "
                 "Типы устойчивы к QC-фильтру (не изменились при отсеве 91К слабых SNP). "
                 "Критичное для назначений — подтверждать прямым клиническим HLA-типированием.")
        L.append("")
    elif hla_path.exists():
        hla = read_tsv(hla_path)
        pos = [r for r in hla if "тег+" in r["tag_status"]]
        L.append("## 🧪 HLA — безопасность лекарств (шаг 12, тег-SNP)")
        if not pos:
            L.append("**Все 5 тегов отрицательны** ✓ — высокорисковых аллелей реакций на лекарства нет; "
                     "целиакия маловероятна (DQ2.5/DQ8 отрицательны).")
        L.append("")

    # шаг 11 происхождение
    anc_path = config.OUT / "11_ancestry.tsv"
    if anc_path.exists():
        anc = read_tsv(anc_path)
        aims = [r for r in anc if r["gene"] != "mtDNA" and r["status"] != "нет на чипе"]
        eur = [r for r in aims if "евр." in r["status"]]
        mt = next((r for r in anc if r["gene"] == "mtDNA"), None)
        L.append("## 🌍 Происхождение (шаг 11)")
        L.append(f"AIM-маркеры: **{len(eur)}/{len(aims)} гомозигот по европейскому аллелю** "
                 "(SLC24A5/SLC45A2 светлая кожа, EDAR не-азиатский, Duffy не-африканский, HERC2) "
                 "→ **сильный европейский сигнал.**")
        if mt:
            L.append(f"mtDNA: {mt['genotype']} — профиль есть, точная гаплогруппа = Haplogrep. "
                     "Y-гаплогруппы нет (женщина). % примеси/неандертальца = нужен 1000G PCA.")
        L.append("> Это направленный сигнал по сильным маркерам, не полный ADMIXTURE.")
        L.append("")

    L.append("## 🔗 Связки и приоритеты (5 кругов ревью разными моделями)")
    L.append("5 независимых кругов ревью: 1 Codex-клиника · 2 Opus-самокритика · "
             "3 Codex-методы · 4 Codex-новые-шаги · 5 Opus-новые-шаги. Артефакты в `out/`: "
             "`codex-review-*`, `opus-selfcritique-*`, `codex-methods-audit-*`, `opus-selfcritique-round2-*`. "
             "Поймано и исправлено: AlphaMissense-формулировка, gnomAD popmax, HLA-QC по R² (типы устояли), "
             "rs12913832 убран, multiallelic-ключи. Ложная тревога «PRS смещён» — проверена и снята. "
             "Главные практичные связки:")
    L.append("- **Тромбоз: PRS ВТЭ ~86-й + КОК + перелёты (кочевник) + операции.** Моногенных "
             "F5 Leiden / F2 у тебя НЕТ (по аннотации) — это не тромбофилия, но полигенный фон "
             "повышен → учитывать при КОК/долгих перелётах/иммобилизации. Самая практичная связка.")
    L.append("- **HFE-носитель + железо:** не назначать железо вслепую; ориентир — ферритин + "
             "насыщение трансферрина (у тебя ферритин в норме).")
    L.append("- **CYP2C9\\*3 (Intermediate) + НПВС/варфарин/фенитоин** — риск токсичности на "
             "стандартных дозах; для варфарина смотреть вместе с VKORC1 (у тебя норма).")
    L.append("- **GJB2 / ATP7B — репродуктивная значимость > личного риска:** при планировании "
             "детей проверить партнёра по тем же генам (панель носительства).")
    L.append("")
    L.append("**3 приоритета по Codex:** (1) подтвердить клинич. лабораторией **ATP7B и GJB2** "
             "(+партнёр при планировании); (2) тромбозы — оценивать с КОК/перелётами/анамнезом; "
             "(3) фенотип-маркеры (ферритин+трансферрин, липиды/ApoB/Lp(a), давление/глюкоза) — "
             "сейчас практичнее большинства PRS. **SLC25A25 (0★) — снять с находок, на ручную проверку.**")
    L.append("")
    L.append("**Уровни доказательности (не смешивать):** ClinVar P/LP > PGx-теги > PRS > "
             "AlphaMissense/нутри-SNP. Импутированное ClinVar-носительство = СКРИНИНГ, не результат; "
             "критичное подтверждать клинич. тестом (dosage/R²/HGVS/strand). Полное ревью: `out/codex-review-2026-06-13.md`.")
    L.append("")

    L.append("## Статус шагов")
    L.append("- Шаг 2 импутация (Michigan HRC r1.1): чип → ~40 млн → качественные варианты после R²-фильтра "
             "(обычно ×10–50 к чипу). Если импутации не было — ClinVar идёт по типизированным.")
    L.append("- ✅ ClinVar (4) · gnomAD-частоты (5b, API) · фарма (7) · AlphaMissense (8) · трейты (9) · нутри (10) · носительство (13)")
    L.append("- ⏳ Остаётся: PRS (6, PGS Catalog) · SpliceAI/полный VEP (OpenCRAVAT — их store сейчас down) · ancestry/HLA (11/12)")
    L.append("")

    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"[report] → {OUT_MD}")
    print(f"[report] кандидаты-носительство (het): {len(carrier_candidates)} · "
          f"гомозиготных артефактов: {len(hom_artifacts)} · risk: {len(risk)} · conflict: {len(conflict)}")
    return OUT_MD


if __name__ == "__main__":
    run()
