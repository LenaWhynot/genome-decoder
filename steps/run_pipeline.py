"""Оркестратор — прогон всего пайплайна по шагам (как pipeline.py у Григория).

  python3 run_pipeline.py

Шаги 1→4→5→7→13→report. Детерминированно; шаг 4 один раз качает ClinVar в .cache/.
"""
import time
import shutil
import subprocess
import config
import step1_normalize
import build_imputed_normalized
import step4_clinvar
import step3_vep
import step5_filter
import step5b_gnomad
import step6_prs
import step7_pharmacogenomics
import step8_alphamissense
import step8b_cadd_spliceai
import step9_traits
import step10_nutrigenomics
import step11_ancestry
import step12_hla
import step13_carrier
import report
from pathlib import Path


def main():
    t0 = time.time()
    print("=" * 64)
    print("ГЕНОМ-ПАЙПЛАЙН — путь: чип → импутация → аннотация")
    print("=" * 64)

    print("\n— ШАГ 1: Нормализация чипа (→ типизированный, по rsID) —")
    step1_normalize.run()   # пишет NORMALIZED_TYPED (для фармы/трейтов/нутри/HLA/ancestry)

    # ШАГ 2: импутация — если есть результаты Michigan, нормализуем их в NORMALIZED (вход
    # для ClinVar/частот/носительства); иначе fallback на типизированный.
    imputed = config.BASE / ".cache" / "imputed" / "chr1.dose.vcf.gz"
    if imputed.exists():
        print("\n— ШАГ 2: Импутация есть → строю импутированный нормализованный —")
        build_imputed_normalized.run()      # пишет NORMALIZED (импутир., chr:pos)
    else:
        print("\n— ШАГ 2: импутации нет → ClinVar по типизированным —")
        shutil.copy(config.NORMALIZED_TYPED, config.NORMALIZED)

    print("\n— ШАГ 4: ClinVar (сверка всех вариантов) —")
    step4_clinvar.run()

    print("\n— ШАГ 5: Фильтр значимости —")
    step5_filter.run()

    print("\n— ШАГ 5b: gnomAD-частоты (сеть) —")
    try:
        step5b_gnomad.run()
    except Exception as e:
        print(f"[step5b] пропуск (нет сети/ошибка API): {e}")

    print("\n— ШАГ 3: VEP-аннотация эффекта (сеть) —")
    try:
        step3_vep.run()
    except Exception as e:
        print(f"[step3] пропуск: {e}")

    print("\n— ШАГ 6: PRS (полигенные риски) —")
    try:
        step6_prs.run()
    except Exception as e:
        print(f"[step6] пропуск (нет сети/ошибка): {e}")

    print("\n— ШАГ 7: Фармакогеномика —")
    step7_pharmacogenomics.run()

    print("\n— ШАГ 8: AlphaMissense (AI-скор) —")
    if (Path(__file__).resolve().parent / ".cache" / "AlphaMissense_hg19.tsv.gz").exists():
        step8_alphamissense.run()
    else:
        print("[step8] пропуск — нет .cache/AlphaMissense_hg19.tsv.gz (скачать 622МБ с Zenodo)")

    print("\n— ШАГ 8b: CADD (+SpliceAI best-effort, сеть) —")
    try:
        step8b_cadd_spliceai.run()
    except Exception as e:
        print(f"[step8b] пропуск: {e}")

    print("\n— ШАГ 9: Трейты —")
    step9_traits.run()

    print("\n— ШАГ 10: Нутригеномика —")
    step10_nutrigenomics.run()

    print("\n— ШАГ 11: Происхождение —")
    step11_ancestry.run()

    print("\n— ШАГ 12: HLA (тег-SNP) —")
    step12_hla.run()

    print("\n— ШАГ 12b: HLA полный (HIBAG, R) —")
    import build_mhc
    build_mhc.run()          # регенерируем MHC-вход с QC по R² (не stale)
    rscript = shutil.which("Rscript") or "/opt/homebrew/bin/Rscript"
    mhc = config.BASE / ".cache" / "mhc_chr6.tsv"
    model = config.BASE / ".cache" / "European-HLA4-hg19.RData"
    if shutil.which(rscript) or (config.BASE / "hla_hibag.R").exists():
        if mhc.exists() and model.exists():
            try:
                subprocess.run([rscript, "hla_hibag.R"], cwd=config.BASE, timeout=600,
                               check=False)
            except Exception as e:
                print(f"[step12b] HIBAG пропуск: {e}")
        else:
            print("[step12b] HIBAG пропуск — нет mhc_chr6.tsv или модели в .cache/")

    print("\n— ШАГ 13: Носительство / ACMG —")
    step13_carrier.run()

    print("\n— Сбор отчёта —")
    report.run()

    print(f"\n✅ Готово за {time.time() - t0:.0f} c. Отчёт в pipeline/out/")


if __name__ == "__main__":
    main()
