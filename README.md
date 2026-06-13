# 🧬 Genome Decoder

> Расшифровка генома за 14 шагов — локально, на чистом Python. Превращает сырой ДНК-тест
> (23andMe / AncestryDNA / Genotek / WGS-VCF) в клинически-осмысленный отчёт: носительство,
> фармакогенетика, полигенные риски, питание под генетику, происхождение, HLA.

**English:** A 14-step local pipeline that decodes a raw genotyping file or WGS VCF into a
human-readable report (carrier status, pharmacogenomics, polygenic risk scores, nutrigenomics,
ancestry, HLA). Pure Python (+ optional R for HLA). **Educational, not medical — confirm any
finding with a clinical lab and a doctor.**

> ⚠️ **Дисклеймер.** Это научпоп-инструмент, **НЕ медицинская диагностика.** Сырые данные
> чипов документированно ложно-завышают редкие патогенные варианты. Любую находку подтверждать
> клинической лабораторией и обсуждать с врачом/генетиком. См. [DISCLAIMER.md](DISCLAIMER.md).

## Что это даёт

- 💊 **Фармакогенетика** — как ваш организм метаболизирует лекарства (CYP2C9/CYP2C19/SLCO1B1/…),
  карточка для врача перед назначением.
- 🧬 **Носительство** — рецессивные болезни (важно при планировании семьи).
- 📊 **Полигенные риски (PRS)** — где вы в популяции по сердцу/диабету/раку/тромбозам.
- 🥗 **Нутригеномика** — лактоза, кофеин, фолат, омега-3, витамины под вашу генетику.
- 🏃 **Трейты** — тип мышц, хронотип, стресс-реактивность.
- 🌍 **Происхождение** + 🧪 **HLA** (drug-safety, целиакия, полное типирование через HIBAG).

## Быстрый старт

```bash
git clone https://github.com/<you>/genome-decoder && cd genome-decoder
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

export GENOME_INPUT="/путь/к/вашему/геному.txt"   # 23andMe/Ancestry .txt или .vcf/.vcf.gz
export SAMPLE_ID="me"
cd steps && python run_pipeline.py
# → отчёт в out/genome-report-<дата>.md
```

Работает с **Claude Code / любым AI-агентом**: положите [SKILL.md](SKILL.md) — агент сам уточнит
формат входа, запустит пайплайн и сделает персональный разбор находок.

## 14 шагов

1. Нормализация (авто-детект 23andMe/VCF) · 2. Импутация (Michigan, опц., ×10–50 покрытия) ·
3. Аннотация эффекта (Ensembl VEP) · 4. ClinVar · 5. gnomAD-фильтр здравого смысла ·
6. GWAS+PRS (PGS Catalog) · 7. Фармакогеномика (CPIC) · 8. AI-скоры (AlphaMissense, CADD) ·
9. Трейты · 10. Нутригеномика · 11. Происхождение (AIM + mtDNA) · 12. HLA (тег-SNP + HIBAG) ·
13. Носительство (ACMG SF) · 14. Сборка отчёта.

Подробно — в [SKILL.md](SKILL.md) (требования, импутация, HLA, ограничения).

## Требования к входу

- **Формат:** 23andMe/AncestryDNA raw `.txt`, или VCF (`.vcf`/`.vcf.gz`), или WGS VCF.
- **Сборка: GRCh37/hg19** (как у raw-чипов). Для GRCh38 нужен лифтовер или замена URL баз.
- Python 3.10+ · ~1–2 ГБ места под базы (ClinVar/AlphaMissense, авто-скачиваются в `.cache/`).

## Принципы (из методологии)

- **Python-детерминизм для механики, AI только для интерпретации.** Каждый вывод — с проверяемым
  источником (ClinVar/gnomAD), не «на глаз».
- **Фильтр здравого смысла:** инделы (чип не типирует) и частые ≥5% (BA1) отсекаются;
  гомозиготные «патогенные» в норме = артефакт референс=минорный-аллель.
- **Уровни доказательности не смешиваются:** ClinVar+gnomAD > PGx/HLA > PRS > AlphaMissense.
- **Мульти-модельное ревью** на надёжность (несколько кругов разными моделями ловят ошибки).

## Credits

Использует открытые ресурсы: ClinVar (NCBI), gnomAD (Broad), Ensembl VEP, AlphaMissense (DeepMind),
CADD, PGS Catalog, Michigan Imputation Server, HIBAG. Спасибо их авторам.

## Лицензия

MIT — см. [LICENSE](LICENSE). Используйте на свой риск; это не медицинский продукт.
