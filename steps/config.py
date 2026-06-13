"""Общие пути и настройки genome-decoder. Чистый stdlib.

Вход задаётся переменными окружения (их выставляет Claude Code по SKILL.md, либо вручную):
  GENOME_INPUT  — путь к файлу генома (23andMe/AncestryDNA .txt ИЛИ VCF; raw или импутированный)
  SAMPLE_ID     — имя образца в выводе (по умолчанию 'SAMPLE')

Авто-детект формата делает step1_normalize. Базы (ClinVar и т.п.) качаются в .cache/.
"""
import os
from pathlib import Path

BASE = Path(__file__).resolve().parent          # .../genome-decoder/steps
ROOT = BASE.parent                              # .../genome-decoder
OUT = ROOT / "out"
CACHE = ROOT / ".cache"
OUT.mkdir(exist_ok=True)
CACHE.mkdir(exist_ok=True)

GENOME_INPUT = Path(os.environ.get("GENOME_INPUT", "")).expanduser()
SAMPLE_ID = os.environ.get("SAMPLE_ID", "SAMPLE")

# нормализованные данные:
# NORMALIZED        — вход для ClinVar/частот/носительства (импутир., если есть импутация)
# NORMALIZED_TYPED  — типизированный по rsID (для фармы/трейтов/нутри/HLA/ancestry)
NORMALIZED = OUT / "01_normalized.tsv"
NORMALIZED_TYPED = OUT / "01_normalized_typed.tsv"

CLINVAR_VCF_GZ = CACHE / "clinvar_GRCh37.vcf.gz"
CLINVAR_HITS = OUT / "04_clinvar_hits.tsv"
SIGNIFICANT = OUT / "05_significant.tsv"
PGX = OUT / "07_pharmacogenomics.tsv"
CARRIER = OUT / "13_carrier.tsv"

# ClinVar GRCh37 (build 37 = совпадает с 23andMe/большинством raw-данных)
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"

# Палиндромные SNP (нельзя однозначно развести по цепи) — помечаем (pre-QC)
AMBIGUOUS_PAIRS = {frozenset("AT"), frozenset("CG")}
