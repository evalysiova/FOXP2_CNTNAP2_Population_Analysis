# Download raw file from repository
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/all-snps.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "all-snps.vcf.gz"))

# Create directories for output plots
output_dir <- file.path("output_snps", "figures")
dir.create(output_dir, recursive = TRUE)

# Required R packages
library("vcfR"); library(adegenet); library(dplyr); library(ggplot2); library(tidyr); library(pegas); library("hierfstat"); library(reshape2)

# Pipeline
vcf <- read.vcfR("data", "raw", "all-snps.vcf.gz")

myID <- getID(vcf)
length(unique(myID, incomparables = NA)) == length(myID)

obj <- vcfR2genind(vcf)
obj

# add population info
url <- "https://github.com/yourusername/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/pop.csv"
download.file(url, destfile = file.path("data", "raw", "pop.csv"))
csv <- read.table(file.path("data", "raw", "pop.csv", header = FALSE, sep = ",")

pop(obj) <- csv[,3]  
obj
objp <- genind2genpop(obj)
objp
popNames(objp)

# Running PCA
pdf(file.path("output_snps", "figures", "PCA_plots_snps.pdf"))
X <- tab(obj, freq=TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

Col1 <- rainbow(length(levels(pop(obj))))
s.class(pca1$li, pop(obj),xax=1,yax=2, col=Col1)
title("PCA of diversity \ axes 1-2")
col <- col <- funky(15)
s.class(pca1$li, fac = pop(obj), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)
dev.off()

# Running CA
pdf(file.path("output_snps", "figures", "CA_plots_snps.pdf"))
ca <- dudi.coa(tab(objp), scannf = FALSE, nf = 3)
ca 
barplot(ca$eig, main = "CA eigenvalues", col = heat.colors(length(ca$eig)))
s.label(ca$li, sub = "CA 1,2", csub = 2)
add.scatter.eig(ca$eig,nf=3,xax=1,yax=3,posi="bottomright")
s.label(ca$li,xax=2,yax=3,lab=popNames(obj),sub="CA 1-3",csub=2)
add.scatter.eig(ca$eig,nf=3,xax=2,yax=3,posi="bottomright")
dev.off()

# Running DAPC
pdf(file.path("output_snps", "figures", "DAPC_plots_snps.pdf"))
# 70 PCs and 5 clusters
grp <- find.clusters(obj, n.pca = 70, n.clust = 5)

table.value(table(pop(obj), grp$grp), col.lab=paste("inf",1:5), row.lab=paste("ori", 1:26))
# 70 PCs and 4 discriminant functions
dapc1 <- dapc(obj, grp$grp, n.pca = 70, n.da = 4)
myCol <- c("orange", "red", "blue", "green", "pink")
scatter(dapc1,xax=1,yax=2, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0, col=Col1,solid=.4, cex=1.5, clab=0, mstree=TRUE, scree.da=FALSE, posi.pca="topright",leg=TRUE,txt.leg=paste("Cluster",1:5))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=2, col=myCol)

myInset <- function()
 {temp <- dapc1$pca.eig; temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black", "lightgray"), c(dapc1$n.pca,1000)), ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance(%)", cex=1, pch=20, type="h", lwd=2) }
  
add.scatter(myInset(), posi="topleft", inset=c(-0.10,-0.10), ratio=.14, bg=transp("white"))
scatter(dapc1, 2, 2, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

compoplot(dapc1,posi="bottomright",txt.leg=paste("Cluster", 1:5), lab=NULL,n.col=1, xlab="individuals")
compoplot(dapc1,subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:5), lab=NULL,n.col=1, xlab="individuals")

myPal <- colorRampPalette(c("blue","gold","red"))

# Identification of variants contributing to population differentiation
# 70 PCs
dapc2 <- dapc(obj, n.pca=70, n.da=2)
scatter(dapc2,xax=1,yax=2, col=transp(myPal(8)), scree.da=FALSE, cell=1.5, cex=1, bg="white",cstar=0, legend=TRUE)
compoplot(dapc2,posi="bottomright",txt.leg=paste("Cluster", 1:5), lab=NULL,n.col=1, xlab="individuals")
compoplot(dapc2,subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:5), lab=NULL,n.col=1, xlab="individuals")
lol <- a.score(dapc2)
par(mar=c(5,4,4,4))
lol <- optim.a.score(dapc2)

# Optimization of DAPC
dapc3 <- dapc(obj, n.pca=50, n.da=2)
temp4 <- which(apply(dapc3$posterior,1, function(e) all(e<0.1)))
lab <- pop(obj)
compoplot(dapc3,subset=temp4, cleg=.5, posi=list(x=0,y=1.2), lab=lab)
x <- dapc3$grp.coord[,1]

group <- rep(1, 5)
group[x > 4.2089248] <- 1
group[x > -1.8099140 & x < -0.6203720] <- 2
group[x > -1.9105798 & x < -1.5057711] <- 3
group[x > -2.3262615 & x < -1.9568285] <- 4
group[x < -3.1228571] <- 5
superpop <- data.frame(x, group)
superpop$group <- c("3", "3", "5", "2", "5", "2", "3", "2", "4", "5", "1", "1", "1", "4", "1", "4", "4", "3", "1", "5", "5", "1", "1", "2", "3", "4")

pal <- colorRampPalette(c("firebrick3", "chartreuse4", "deeppink3", "cadetblue2", "darkgoldenrod2"))
palette(pal(5))
scatter(dapc3,xax=1,yax=2, col= superpop$group, solid=.9, scree.da=FALSE,cell=1.5, cex=1.5,  bg="white",cstar=0, legend=TRUE)
scatter(dapc3,xax=1,yax=2, col= superpop$group, solid=.9, scree.da=FALSE,cell=2.5, cex=1.5, bg="white",cstar=0, legend=TRUE)
set.seed(4)

contrib <- loadingplot(dapc3$var.contr, axis=2, thres=0.0008, lab.jitter=1)
dev.off()

contrib_snps_sorted <- readLines("contrib_snps_sorted.txt")
pdf(file.path("output_snps", "figures", "all_freq_snps.pdf"), width = 14, height = 6)
for (id in contrib_snps_sorted) 
 {snp <- tab(genind2genpop(obj[loc = c(id)]), freq=TRUE)
  all_alleles <- data.frame(ID = getID(vcf), REF = getREF(vcf), ALT = getALT(vcf))
  if (id %in% all_alleles$ID)
    {matching <- all_alleles[all_alleles$ID == id, ]
      print(matching)
      alleles <- c(matching$REF, matching$ALT)}
      matplot(snp, pch=1:length(alleles), type="b", xlab="Population", ylab="Allele Frequency", xaxt="n",cex=1.5, main=paste(id))
          axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
          legend("bottomright", legend = alleles, pch = 1:length(alleles), pt.cex = 1.5, bty = "n")
  }
dev.off()

#F-statistics and Heterozygosity
fstab <- Fst(as.loci(obj))
Hs(obj)

Hobs <- lapply(seppop(obj), function(e) mean(summary(e)$Hobs, na.rm = TRUE))
populationn <- c("GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV","ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT","LWK", "ASW", "MXL", "TSI", "GIH")
Observed <- c()
Expected <- c()
df <- data.frame(populationn,Observed,Expected)
print(df)
combined <- rbind(Observed = df$Observed, Expected = df$Expected)
obj_sum <- summary(obj)
Ho<-mean(obj_sum$Hobs)
He<-mean(obj_sum$Hexp)
F <- (He - Ho)/He
result <- (df$Expected - df$Observed) / df$Expected
size <- table(pop(obj))
sample_size <- c(91,99,105,104,93,94,107,85,96,99,96,113,99,86,85,102,102,99,108,103,104,99,61,64,107,103)
populations <- seppop(obj)
alleles_per_pop <- sapply(populations, function(e) sum(nAll(e)))

F_value <- c()

df1 <- data.frame(populationn,Observed,Expected,sample_size,alleles_per_pop,F_value)
write.table(df1, file = file.path("output_snps", "~/dataframe2.tsv"), sep = \"t")

contrib_snps_ids <- c("rs2245192", "rs7782412", "rs7799652", "rs7799269", "rs72611587", "rs1603450", "rs58525662", "rs10250389", "rs1496545", "rs6962971", "rs885610", "rs12670824", "rs6975114", "rs60424529", "rs17170378", "rs12538316", "rs17170379", "rs10276382", "rs10279700", "rs28840228", "rs6974852", "rs12534978", "rs1525215", "rs10261539", "rs10248220", 
"rs10085766", "rs10256646", "rs1383011", "rs1021008", "rs2692182", "rs899613", "rs899615", "rs2692176", "rs1917600", "rs2692157", "rs10254102", "rs2620454", "rs201667516", "rs10485845", "rs2707559", "rs10485844", "rs7456839", "rs10243597", "rs4726920", "rs12671444", "rs6961656", "rs10441210", "rs7797715", "rs2040920", "rs4726922", "rs6464853", "rs12669189",
"rs6464854", "rs7791252", "rs6962625", "rs2214722", "rs6964783", "rs4725763", "rs4726923", "rs17170838", "rs7803754", "rs7781615", "rs1859540", "rs1859539")

contrib_obj <- obj[loc = contrib_snps_ids]

statistics <- basic.stats(contrib_obj)
Ho_total <- statistics$Ho
Hexp_total <- statistics$Hs
mean_Ho <- as.numeric(rowMeans(Ho_total, na.rm = TRUE))
mean_Hexp <- as.numeric(rowMeans(Hexp_total, na.rm = TRUE))

all_allele_freq <- tab(contrib_obj, freq = TRUE)
total_freq <- colMeans(all_allele_freq, na.rm = TRUE)
freq_df <- data.frame(Allele = names(total_freq), Frequency = as.numeric(total_freq), stringsAsFactors = FALSE)

total_ref_freq <- c()
total_alt_freq <- c()
for(i in 1:nrow(freq_df)){
if(grepl("\\.0$", freq_df$Allele[i]))
        {total_ref_freq <- c(total_ref_freq, freq_df$Frequency[i])}
else if(grepl("\\.1$", freq_df$Allele[i]))
        {total_alt_freq <- c(total_alt_freq, freq_df$Frequency[i])}}
        
# Download raw file (contributing SNPs) from repository
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/contrib_snps.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "contrib_snps.vcf.gz"))

contrib_vcf <- read.vcfR(file.path("data", "raw", "contrib_snps.vcf.gz")
ref_alleles <- getREF(contrib_vcf)
alt_all <- getALT(contrib_vcf)

alt_alleles <- c()
for(i in 1:length(alt_all))
        {if(grepl("\\,", alt_all[i]))
                {alt_alleles <- c(alt_alleles, gsub(",.*", "", alt_all[i]))}
        else
                {alt_alleles <- c(alt_alleles, alt_all[i]) }}



total_statistics_df <- data.frame(IDS = contrib_snps_ids, Ref = ref_alleles, Alt = alt_alleles, Ref_frequency = total_ref_freq, Alt_freq = total_alt_freq, Ho = mean_Ho, He = mean_Hexp)
write.table(total_statistics_df, file = file.path("output_snps", "total_per_locus2.tsv"), sep = \"t")

stat_per_pop_df <- data.frame(ID = contrib_snps_ids)
for(i in 1:length(populationn))
        {name <- populationn[i]
        variable_name1 <- paste0("Ho_", name); variable_name2 <- paste0("He_", name)
        assign(variable_name1, as.numeric(Ho_total[, name]))
        assign(variable_name2, as.numeric(Hexp_total[, name]))
        stat_per_pop_df[[variable_name1]] <- get(variable_name1)
        stat_per_pop_df[[variable_name2]] <- get(variable_name2)}

write.table(stat_per_pop_df, file = file.path("output_snps", "stat_per_pop2.tsv"), sep = "\t")

long_het <- stat_per_pop_df %>% pivot_longer(-ID, names_to = c(".value", "Population"), names_pattern = "(Ho|He)_(.+)")
fis_df <- long_het %>% mutate(Fis = 1-(Ho/He))
shapiro.test(fis_df$Fis)
# 
# Shapiro-Wilk normality test
# 
# data:  fis_df$Fis
# W = 0.98463, p-value = 2.413e-12

fis_results <- fis_df %>% group_by(Population) %>% summarise(mean_Fis = mean(Fis, na.rm = TRUE), p_value = wilcox.test(na.omit(Fis), mu = 0, exact = FALSE)$p.value, 
		 W_stat = wilcox.test(na.omit(Fis), mu = 0)$statistic) %>%
  mutate(direction = ifelse(mean_Fis > 0, "Inbreeding", "Excess Heterozygosity"), sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value <0.05 ~ "*"), 
         label_pos = ifelse(mean_Fis >= 0, mean_Fis + 0.003, mean_Fis - 0.003)) %>%
  arrange(p_value)
  
pdf(file.path("output_snps", "figures", "Fis_of_contributing_snps.pdf"))
ggplot(fis_results, aes(x = reorder(Population, mean_Fis), y = mean_Fis, fill = direction)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(y = label_pos, label = sig), hjust = 0.5, size  = 3.5, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Inbreeding" = "salmon", "Excess Heterozygosity" = "turquoise3")) +
  labs(title = "Mean Inbreeding Coefficient (Fis) per Population", subtitle = "Wilcoxon test | * p<0.05  ** p<0.01 *** p<0.001", x = "Population", y = "Mean Fis", fill = "direction") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey40", size = 8))
dev.off()

long_df <- stat_per_pop_df %>% pivot_longer(-ID, names_to = c(".value", "Population"), names_pattern = "(Ho|He)_(.+)")

superpop <- data.frame(Population = c("ACB","ASW","ESN","GWD","LWK","MSL","YRI", "CEU","FIN","GBR","IBS","TSI", "CDX","CHB","CHS","JPT","KHV", "BEB","GIH","ITU","PJL","STU","CLM","MXL","PEL","PUR"), Superpop = c(rep("AFR", 7), rep("EUR", 5), rep("EAS", 5), rep("SAS", 5), rep("AMR", 4)))

pop_order <- superpop$Population

pdf(file.path("output_snps", "figures", "Ho-He_of contributing_snps.pdf"))
long_df %>% mutate(diff = Ho - He) %>% left_join(superpop, by = "Population") %>%
  mutate(Population = factor(Population, levels = pop_order)) %>%
  ggplot(aes(x = Population, y = diff, color = ID)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 0.8, linetype  = "solid") +
  annotate("rect", xmin = 0.5, xmax = 7.5, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = 0.07) +
  annotate("rect", xmin = 7.5, xmax = 12.5, ymin = -Inf, ymax = Inf, fill = "#F0E442", alpha = 0.07) +
  annotate("rect", xmin = 12.5, xmax = 17.5, ymin = -Inf, ymax = Inf, fill = "#009E73", alpha = 0.07) +
  annotate("rect", xmin = 17.5, xmax = 22.5, ymin = -Inf, ymax = Inf, fill = "#CC79A7", alpha = 0.07) +
  annotate("rect", xmin = 22.5, xmax = 26.5, ymin = -Inf, ymax = Inf, fill = "#56B4E9", alpha = 0.07) +
  annotate("text", x = 4, y = Inf, label = "AFR", vjust = 2, fontface = "bold", size = 3.5, color = "#E69F00") +
  annotate("text", x = 10, y = Inf, label = "EUR", vjust = 2, fontface = "bold", size = 3.5, color = "#F0E442") +
  annotate("text", x = 15, y = Inf, label = "EAS", vjust = 2, fontface = "bold", size = 3.5, color = "#009E73") +
  annotate("text", x = 20, y = Inf, label = "SAS", vjust = 2, fontface = "bold", size = 3.5, color = "#CC79A7") +
  annotate("text", x = 24.5, y = Inf, label = "AMR", vjust = 2, fontface = "bold", size = 3.5, color = "#56B4E9") +
  labs(title = "Ho - He per INDEL per Population", x = "Population", y = "Ho - He", color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey40"))
dev.off()

# Pairwise Fst
pairwise_fst <- pairwise.neifst(genind2hierfstat(contrib_obj))
fst_matrix <- as.matrix(pairwise_fst)
superpop_map <- c(ACB="AFR", ASW="AFR", ESN="AFR", GWD="AFR", LWK="AFR", MSL="AFR", YRI="AFR",
  CEU="EUR", FIN="EUR", GBR="EUR", IBS="EUR", TSI="EUR",CDX="EAS", CHB="EAS", CHS="EAS", JPT="EAS",
  KHV="EAS", CLM="AMR", MXL="AMR", PEL="AMR", PUR="AMR", BEB="SAS", GIH="SAS", ITU="SAS", PJL="SAS",
  STU="SAS")

pop <- rownames(fst_matrix); superpop <- superpop_map[pop]
order_i <- pop[order(superpop)]
fst_matrix <- fst_matrix[order_i, order_i]
fst_long <- reshape2::melt(fst_matrix, na.rm = FALSE); colnames(fst_long) <- c("Pop1", "Pop2", "Fst")
fst_long$Pop1 <- factor(fst_long$Pop1, levels = order_i)
fst_long$Pop2 <- factor(fst_long$Pop2, levels = order_i)
superpop_df <- data.frame(pop = factor(order_i, levels = order_i), superpop = superpop_map[order_i])
colors_superp <- c("AFR" = "#E69F00", "AMR" = "#56B4E9", "EUR" = "#F0E442", "EAS" = "#009E73", "SAS" = "#CC79A7")

pdf(file.path("output_snps", "figures", "Pairwise_Fst_of_contributing_snps.pdf"))
heatmap <- ggplot(fst_long, aes(x = Pop1, y = Pop2, fill = Fst)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightpink", high = "darkred", limits = c(-0.01, 0.6), na.value = "grey90", name = expression(F[ST])) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7), panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.05), legend.position = "right")
guide_left <- ggplot(superpop_df, aes(x = 1, y = pop, fill = superpop)) +
  geom_tile(color = "white")+
  scale_fill_manual(values = colors_superp, name = "Superpopulation") + 
  theme_void() +
  theme(legend.position = "none")
guide_bottom <- ggplot(superpop_df, aes(x = pop, y = 1, fill = superpop)) +
  geom_tile(color = "white")+
  scale_fill_manual(values = colors_superp, name = "Superpopulation") + 
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 8))
(guide_left + heatmap + plot_layout(widths = c(0.03, 1))) / 
  (plot_spacer() + guide_bottom + plot_layout(widths = c(0.03, 1))) +
  plot_layout(heights = c(1, 0.03))
dev.off()


# Mean Fst - Ho correlation
mean_fst <- function(fst_matrix)
{diag(fst_matrix )<- NA
data.frame(Population = rownames(fst_matrix), mean_Fst = rowMeans(fst_matrix, na.rm = TRUE))
}

fst_long <- mean_fst(pairwise_fst)

Ho_total <- as.data.frame(Hobs)
Ho_total$SNP <- rownames(Ho_total)
ho_long <- Ho_total %>% pivot_longer(cols = -SNP, names_to = "Population", values_to = "Ho")
correlation <- left_join(ho_long, fst_long, by = c("Population"))
correl_results <- correlation |> summarise(rho = cor(Ho, mean_Fst, method = "spearman", use = "complete.obs"), p_val = cor.test(Ho, mean_Fst, method = "spearman", exact = FALSE)$p.value, .groups = "drop") |> mutate (significant = ifelse(p_val <- 0.05, "*", "ns"))

pdf(file.path("output_snps", "figures", "Fst_Ho_correlation_snps.pdf"))
ggplot(correlation, aes(x = Ho, y = mean_Fst, color = Superpopulation, label = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed", linewidth = 0.7) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  scale_color_manual(values = c("AFR" = "#E69D00", "AMR" = "#56B4E9", "EUR" = "#F0E442", "EAS" = "#009E73", "SAS" = "#CC79A7")) + 
  labs(title = expression("Correlation Between Ho and mean Fst"),
       x = expression("Observed Heterozygosity (Ho)"), y = expression("Mean pairwise Fst")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position  = "bottom", strip.text = element_text(face = "bold"))
dev.off()

