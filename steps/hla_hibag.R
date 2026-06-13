#!/usr/bin/env Rscript
# ШАГ 12 (полное HLA) — HIBAG-типирование из импутированных SNP.
# Модель: European-HLA4-hg19 (A/B/C/DRB1/DQA1/DQB1/DPB1). Матч по позиции (импут. без rsID).
# Вывод: out/12_hibag.tsv (HLA-аллели + вероятность).
suppressMessages(library(HIBAG))

load(".cache/European-HLA4-hg19.RData")          # объект HLA4
geno <- read.table(".cache/mhc_chr6.tsv", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
geno <- geno[!duplicated(geno$pos), ]            # одна запись на позицию
cat(sprintf("[hibag] MHC SNP (импутир.): %d\n", nrow(geno)))

# genoObj: genotype = число копий ПЕРВОГО аллеля (ref) в snp.allele "ref/alt";
# наш geno = число копий ALT → ref-копии = 2 - geno. snp.id = позиция (для матча).
snp.id     <- as.character(geno$pos)
snp.allele <- paste(geno$ref, geno$alt, sep = "/")
genotype   <- matrix(as.integer(2L - geno$geno), ncol = 1)

genoObj <- list(genotype = genotype, sample.id = Sys.getenv("SAMPLE_ID", "SAMPLE"),
             snp.id = snp.id, snp.position = geno$pos,
             snp.allele = snp.allele, assembly = "hg19")
class(genoObj) <- "hlaSNPGenoClass"

res <- data.frame()
for (gene in names(HLA4)) {
  mobj <- HLA4[[gene]]
  mobj$snp.id <- as.character(mobj$snp.position)   # переименовываем под позиционный матч
  ov <- sum(mobj$snp.id %in% snp.id)
  model <- tryCatch(hlaModelFromObj(mobj), error = function(e) NULL)
  if (is.null(model)) next
  pred <- tryCatch(predict(model, genoObj, type = "response", match.type = "Position"),
                   error = function(e) { cat(sprintf("[hibag] %s: %s\n", gene, conditionMessage(e))); NULL })
  if (is.null(pred)) next
  v <- pred$value
  a1 <- v$allele1[1]; a2 <- v$allele2[1]; pr <- v$prob[1]
  cat(sprintf("[hibag] HLA-%-5s %s / %s   (prob %.2f, SNP-overlap %d/%d)\n",
              gene, a1, a2, pr, ov, length(mobj$snp.id)))
  res <- rbind(res, data.frame(gene = paste0("HLA-", gene),
               allele1 = a1, allele2 = a2, prob = round(pr, 3),
               snp_overlap = ov, snp_model = length(mobj$snp.id)))
}

write.table(res, "out/12_hibag.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("[hibag] → out/12_hibag.tsv\n")
