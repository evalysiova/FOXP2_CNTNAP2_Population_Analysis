# Download raw file from repository
dir.create(file.path("data", "raw"), recursive = TRUE)
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/common_merged.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "common_merged.vcf.gz"))

output_dir <- file.path("output_LD", "figures")
dir.create(output_dir, recursive = TRUE)

# Required R packages
library(vcfR); library(GWLD); library(dplyr); library(ggplot2)

vcf <- read.vcfR(file.path("data", "raw", "common_merged.vcf.gz"))
genotypes <- extract.gt(vcf, element = "GT")

convert_gt <- function(x) 
{x <- gsub("\\|", "/", x)
 sapply(x, function(g) 
 {if(is.na(g) || g == "./." || g == ".|.") return(NA)
  alleles <- strsplit(g, "/")[[1]]
  alleles_num <- sort(as.numeric(alleles))
  sum(as.numeric(alleles))
  if(all(alleles_num <= 1))
    return(sum(alleles_num))
  paste <- paste(alleles_num, collapse = "/")
  switch(paste, "0/2" = 3, "1/2" = 4, "2/2" = 5)
})
}

geno_num <- apply(genotypes, 2, convert_gt)
row.names(geno_num) <- row.names(genotypes)
colnames(geno_num) <- colnames(genotypes)
geno_num <- t(geno_num)

pop_clean <- read.table("pop.csv", header = FALSE, sep = ",")
Info <- data.frame(CHROM = getCHROM(vcf), POS   = getPOS(vcf), ID    = getID(vcf))
mycols <- colorRampPalette(c("#ffe6e6", "#ff0000", "#990000"))(100)

superpop_map <- c(
  ACB="AFR", ASW="AFR", ESN="AFR", GWD="AFR", LWK="AFR", MSL="AFR", YRI="AFR",
  CEU="EUR", FIN="EUR", GBR="EUR", IBS="EUR", TSI="EUR",
  CDX="EAS", CHB="EAS", CHS="EAS", JPT="EAS", KHV="EAS",
  CLM="AMR", MXL="AMR", PEL="AMR", PUR="AMR",
  BEB="SAS", GIH="SAS", ITU="SAS", PJL="SAS", STU="SAS")

pop_clean <- pop_clean %>% rename(id = V1, sample = V2, pop = V3) %>% mutate(superpop = superpop_map[pop])

superpop_unique <- c("AFR", "EUR", "EAS", "SAS", "AMR")
rmi_list <- list()
target_snp <- "rs7782412"


pdf(file.path("output_LD", "figures", "LD_trait.pdf"), width = 9,  height = 6)
for(i in seq_along(superpop_unique))
{spop_i <- superpop_unique[[i]]
ids <- pop_clean[pop_clean[, 4] ==spop_i, ]
geno_num_i <- geno_num[row.names(geno_num) %in% ids[, 2], ]
snp_labels <- colnames(geno_num_i)
snp_sizes <- ifelse(snp_labels %in% target_snp, 10, 0.001)
rmi <- GWLD(geno_num_i, method = "RMI", cores = 4)
rmi_list[[spop_i]] <- GWLD(geno_num_i, method = "RMI", cores = 4)
p <- HeatMap(geno_num_i, method = "RMI", SnpPosition = Info[, 1:2],
  			SnpName = snp_labels, cores = 4, color = "blueTored", label.size = snp_sizes)

grid::grid.text(spop_i, x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 14, fontface = "bold", col = "black"))
}

dev.off()

pdf(file.path("output_LD", "figures", "LD_trait_heatmap.pdf"))
for(i in seq_along(rmi_list))
  {rmi_full <- rmi_list[[i]]
   rmi_full[lower.tri(rmi_full)] <- t(rmi_full)[lower.tri(rmi_full)]
   colnames(rmi_full) <- as.character(Info$ID)
   rownames(rmi_full) <- as.character(Info$ID)
   rmi_full2 <- rmi_full[31:89, 31:89]
   rmi_full_sub <- rmi_full[51:69, 51:69]
   rmi_long <- reshape2::melt(rmi_full2, na.rm = FALSE)
   rmi_long2 <- reshape2::melt(rmi_full_sub, na.rm = FALSE)
   colnames(rmi_long2) <- c("VAR1", "VAR2", "RMI"); colnames(rmi_long) <- c("VAR1", "VAR2", "RMI")
   p <- ggplot(rmi_long, aes(x = VAR1, y = VAR2, fill = RMI)) +
     geom_tile(color = "white") +
     scale_fill_gradient(low = "#FFE4E1", high = "darkred", limits = c(-0.1, 1.5), na.value = "gray") +
     coord_fixed() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
   q <- ggplot(rmi_long2, aes(x = VAR1, y = VAR2, fill = RMI)) +
     geom_tile(color = "white") +
     scale_fill_gradient(low = "#FFE4E1", high = "darkred", limits = c(-0.1, 1.5), na.value = "gray") +
     coord_fixed() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
   print(p); print(q)}
dev.off()

