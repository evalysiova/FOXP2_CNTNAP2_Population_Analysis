# Download raw file from repository
dir.create(file.path("data", "raw"), recursive = TRUE)
url <- "https://github.com/yourusername/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/all-indels.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "all-indels.vcf.gz"))

# Create directories for output plots
output_dir <- file.path("output_indels", "figures")
dor.create(output_dir, recursive = TRUE)

# Required R packages
library("vcfR"); library(adegenet); library(ggplot2); library(reshape2); library(dplyr); library(pegas); library(hierfstat); library(tidyr)

# Pipeline
vcf <- read.vcfR("data/raw/all-indels.vcf.gz")

myID <- getID(vcf)
length(unique(myID, incomparables = NA)) == length(myID)
duplicates <- which(duplicated(myID) | duplicated(myID, fromLast = TRUE))
dup_ids <- vcf@fix[duplicates, c("CHROM","POS","ID")]

fix <- vcf@fix
fix[,"ID"] <- make.unique(fix[,"ID"], sep = "_")
vcf@fix <- fix
# check again for duplicates
myID2 <- getID(vcf)
length(unique(myID2, incomparables = NA)) == length(myID2)

obj <- vcfR2genind(vcf)

# add population info
url <- "https://github.com/yourusername/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/pop.csv"
download.file(url, destfile = file.path("data", "raw", "pop.csv"))
csv <- read.table(file.path("data", "raw", "pop.csv", header = FALSE, sep = ",")
# check order of IDs
inds <- indNames(obj)
all(inds == csv$V3)

pop(obj) <- csv[,3]  
obj
objp <- genind2genpop(obj)
objp
popNames(objp)

# Running PCA
pdf(file.path("output_indels", "figures", "PCA_plots_indels.pdf"))
X <- tab(obj, freq=TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

Col1 <- rainbow(length(levels(pop(obj))))
s.class(pca1$li, pop(obj),xax=1,yax=2, col=Col1)
title("PCA of diversity \ axes 1-2")
col <- funky(15)
s.class(pca1$li, fac = pop(obj), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)
colorplot(pca1$li[c(1,2)], pca1$li, transp=TRUE, cex=1.5, xlab="PC 1", ylab="PC 2")
title("PCA of human genes \ axes 1-2")
abline(v=0,h=0,col="grey", lty=2)
dev.off()

# Running CA
pdf(file.path("output_indels", "figures", "CA_plots_indels.pdf"))
ca <- dudi.coa(tab(objp), scannf = FALSE, nf = 3)
ca 
barplot(ca$eig, main = "CA eigenvalues", col = heat.colors(length(ca$eig)))
s.label(ca$li, sub = "CA 1,2", csub = 2)
add.scatter.eig(ca$eig,nf=3,xax=1,yax=3,posi="bottomright")
s.label(ca$li,xax=2,yax=3,lab=popNames(obj),sub="CA 1-3",csub=2)
add.scatter.eig(ca$eig,nf=3,xax=2,yax=3,posi="bottomright")
dev.off()

# Running DAPC
pdf(file.path("output_indels", "figures", "DAPC_plots_indels.pdf"))
# 100 PCs and 5 clusters
grp <- find.clusters(obj, n.pca = 100, n.clust = 5)

table.value(table(pop(obj), grp$grp), col.lab=paste("inf",1:5), row.lab=paste("ori", 1:26))
# 100 PCs and 4 discriminant functions
dapc1 <- dapc(obj, grp$grp, n.pca = 100, n.da = 4)
myCol <- c("orange", "red", "blue", "green", "pink")
scatter(dapc1,xax=1,yax=2, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0, col=Col1,solid=.4, cex=1.5, clab=0, mstree=TRUE, scree.da=FALSE, posi.pca="topright",leg=TRUE,txt.leg=paste("Cluster",1:5))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=2, col=myCol)

myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black", "lightgray"),
                     c(dapc1$n.pca,1000)), ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance
(%)",
       cex=1, pch=20, type="h", lwd=2) }
add.scatter(myInset(), posi="topleft", inset=c(-0.10,-0.10), ratio=.14, bg=transp("white"))
scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)
myPal <- colorRampPalette(c("blue","gold","red"))

# Identification of variants contributing to population differentiation
# 60 PCs
dapc2 <- dapc(obj, n.pca=60, n.da=2)
scatter(dapc2,xax=1,yax=2, col=transp(myPal(8)), scree.da=FALSE, cell=1.5, cex=1, bg="white",cstar=0, legend=TRUE)
compoplot(dapc1,posi="bottomright",txt.leg=paste("Cluster", 1:5), lab=NULL,n.col=1, xlab="individuals")
compoplot(dapc1,subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:5), lab=NULL,n.col=1, xlab="individuals")

# Optimization of DAPC
lol <- a.score(dapc2)
par(mar=c(5,4,4,4))
lol <- optim.a.score(dapc2)
dapc3 <- dapc(obj, n.pca=60, n.da=2)
temp4 <- which(apply(dapc3$posterior,1, function(e) all(e<0.1)))
lab <- pop(obj)
compoplot(dapc2,subset=temp4, cleg=.5, posi=list(x=0,y=1.2), lab=lab)
x <- dapc3$grp.coord[,1]
group <- rep(1, 5)
group[x > 4.3137145] <- 1
group[x > -1.7658436 & x < -0.5433659] <- 2
group[x > -1.8969794 & x < -1.4355424] <- 3
group[x > -2.2897167 & x < -1.8400496] <- 4
group[x < -3.4463822] <- 5
superpop <- data.frame(x, group)
superpop$group <- c("3", "3", "5", "2", "5", "2", "3", "2", "4", "5", "1", "1", "1", "4", "1", "4", "4", "3", "1", "5", "5", "1", "1", "2", "3", "4")

pal <- colorRampPalette(c("firebrick3", "chartreuse4", "deeppink3", "cadetblue2", "darkgoldenrod2"))
palette(pal(5))
scatter(dapc3,xax=1,yax=2, col= superpop$group, solid=.9, scree.da=FALSE,cell=1.5, cex=1.5,
        +         bg="white",cstar=0, legend=TRUE)
scatter(dapc3,xax=1,yax=2, col= superpop$group, solid=.9, scree.da=FALSE,cell=2.5, cex=1.5,
        +         bg="white",cstar=0, legend=TRUE)


set.seed(4)
contrib <- loadingplot(dapc3$var.contr, axis=2, thres=0.004, lab.jitter=1, lab=" ")
dev.off()

pdf(file.path("output_indels", "figures", "all_freq_indels.pdf"), width = 14, height = 6)
#INDEL 1 => rs10585289
indel1 <- tab(genind2genpop(obj[loc=c("rs10585289")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel1, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs10585289")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("GTTTA", "G"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL2 => rs34526725
indel2 <- tab(genind2genpop(obj[loc=c("rs34526725")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel2, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs34526725")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("GA", "G"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL3 => rs71999274 !3 ALLELES!
indel3 <- tab(genind2genpop(obj[loc=c("rs71999274")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel3, pch=c(1,2,3), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs71999274 & rs7799652")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("topright", legend = c("TA", "G", "GA"), pch = c(1, 2, 3), pt.cex = 1.5, bty = "n")
matplot(indel3[, 2:3], pch=c(2,3), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs71999274")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("topright", legend = c("G", "GA"), pch = c(2, 3), pt.cex = 1.5, bty = "n")

#INDEL4 => rs143319023 !3 alleles!
indel4 <- tab(genind2genpop(obj[loc=c("rs143319023")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel4, pch=c(1,2,3), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs143319023 ")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("topright", legend = c("AT", "ATT", "A"), pch = c(1, 2, 3), pt.cex = 1, bty = "n")

#INDEL5 => rs10553315
indel5 <- tab(genind2genpop(obj[loc=c("rs10553315")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel5, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs10553315")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("TAA", "T"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDELS6 => rs34910578
indel6 <- tab(genind2genpop(obj[loc=c("rs34910578")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel6, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs34910578")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("TAAA", "T"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL7 => rs202106603 !3 ALLELES!
indel7 <- tab(genind2genpop(obj[loc=c("rs202106603")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel7, pch=c(1,2,3), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs202106603")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("topright", legend = c("AT", "A", "TT"), pch = c(1, 2, 3), pt.cex = 1.5, bty = "n")

#INDEL8 => rs59499448
indel8 <- tab(genind2genpop(obj[loc=c("rs59499448")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel8, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs59499448")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("CCT", "C"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL9 => rs542276531
indel9 <- tab(genind2genpop(obj[loc=c("rs542276531")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel9, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs542276531")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("AAT", "A"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL10 => rs556222360
indel10 <- tab(genind2genpop(obj[loc=c("rs556222360")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel10, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs556222360")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("A", "AC"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL11 => rs200974074
indel11 <- tab(genind2genpop(obj[loc=c("rs200974074")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel11, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs200974074")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("CAA", "C"), pch = c(1, 2), pt.cex = 1.5, bty = "n")

#INDEL12 => rs58192570
indel12 <- tab(genind2genpop(obj[loc=c("rs58192570")]),freq=TRUE)
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,3.0),las=3)
matplot(indel12, pch=c(1,2), type="b",xlab="population",ylab="allele frequency", xaxt="n",cex=1.5, main="rs58192570")
axis(side = 1, at = 1:26, labels = popNames(obj), cex.axis = 0.6)
legend("bottomright", legend = c("G", "GA"), pch = c(1, 2), pt.cex = 1.5, bty = "n")
dev.off()

#F-statistics and Heterozygosity
fstab <- Fst(as.loci(obj))
Hs(obj)

Hobs <- lapply(seppop(obj), function(e) mean(summary(e)$Hobs, na.rm = TRUE))

Observed <- c(0.1041878, 0.1016299,0.08859358,0.1159854, 0.08861251,0.1129466,0.1052024,0.1006675,0.1032852,0.09086381,0.1343368,0.1355775,0.1353433,0.1063193,0.1380482,0.1058955,0.104316,0.1028299,0.1345433,0.08819473,0.09158459,0.1330594,0.1366549,0.1072823,0.10504,0.103782)
Expected <- c(0.10263148,0.10048418,0.08838874,0.11230178,0.08923160,0.11078311,0.10222575,0.09879175,0.10537949,0.08951976,0.13235790,0.13327245,0.13055418,0.10228623,0.13393897,0.10311290,0.10168254,0.10251496,0.13254440,0.08808437,0.08955418,0.13107118,0.13391530,0.10689850,0.10266974,0.10200252)
df <- data.frame(populationn,Observed,Expected)
print(df)
combined <- rbind(Observed = df$Observed, Expected = df$Expected)
pdf("plott.pdf", width = 17, height = 8)
barplot(combined, beside = TRUE, col = c("skyblue", "salmon"), legend.text = TRUE, args.legend = list(x = "topright", cex = 0.8), names.arg = df$populationn, main = "Observed vs Expected Heterozygosity", xlab = "Population", ylab = "Heterozygosity")
obj_sum <- summary(obj)
Ho<-mean(obj_sum$Hobs)
He<-mean(obj_sum$Hexp)
F <- (He - Ho)/He
result <- (df$Expected - df$Observed) / df$Expected
size <- table(pop(obj))
sample_size <- c(91,99,105,104,93,94,107,85,96,99,96,113,99,86,85,102,102,99,108,103,104,99,61,64,107,103)
populations <- seppop(obj)
alleles_per_pop <- sapply(populations, function(e) sum(nAll(e)))
F_value <- c(-0.015164158,-0.011401994,-0.002317490,-0.032801083,0.006938013,-0.019529060,-0.029118397,-0.018986909,0.019873791,-0.015014004,-0.014951129,-0.017295773,-0.036683008,-0.039429257,-0.030679869,-0.026985954,-0.025898842,-0.003072137,-0.015080984,-0.001252890,-0.022672420,-0.015169010,-0.020457707,-0.003590322,-0.023086257,-0.017445451)
df1 <- data.frame(populationn,Observed,Expected,sample_size,alleles_per_pop,F_value)
write.csv(df1, file = "~/dataframe.xlsl")

contrib_indel_ids <- c("rs10585289", "rs34526725", "rs71999274", "rs143319023", "rs10553315", "rs34910578", "rs202106603", "rs59499448", "rs542276531", "rs556222360", "rs200974074", "rs58192570")
#populationn already defined

contrib_obj <- obj[loc = contrib_indel_ids]

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

total_alt2_freq <- c("-", "-", "0.049836928", "0.003194888", "-", "-", "0.008935703", "-", "-", "-", "-", "-")

# Download raw file (contributing SNPs) from repository
url <- "https://github.com/yourusername/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/contrib_indels.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "contrib_indels.vcf.gz"))

contrib_vcf <- read.vcfR(file.path("data", "raw", "contrib_indels.vcf.gz")
ref_alleles <- getREF(contrib_vcf)
alt_all <- getALT(contrib_vcf)

alt_alleles <- c()
for(i in 1:length(alt_all))
        {if(grepl("\\,", alt_all[i]))
                {alt_alleles <- c(alt_alleles, gsub(",.*", "", alt_all[i]))}
        else
                {alt_alleles <- c(alt_alleles, alt_all[i]) }}

alt2_alleles <- c("-", "-", "G", "ATT", "-", "-", "A", "-", "-", "-", "-", "-")


total_statistics_per_locus_df <- data.frame(ID = contrib_indel_ids, Ref = ref_alleles, Alt1 = alt_alleles, Alt2 = alt2_alleles, Ref_Frequency = total_ref_freq, Alt1_Frequency = total_alt_freq, Alt2_Frequency = total_alt2_freq, Ho = mean_Ho, He = mean_Hexp)
write.table(total_statistics_per_locus_df, file = file.path("output_indels", "total_per_locus.tsv"), sep = "\t")

stat_per_pop_df <- data.frame(ID = contrib_indel_ids)
for(i in 1:length(populationn))
        {name <- populationn[i]
        variable_name1 <- paste0("Ho_", name); variable_name2 <- paste0("He_", name)
        assign(variable_name1, as.numeric(Ho_total[, name]))
        assign(variable_name2, as.numeric(Hexp_total[, name]))
        stat_per_pop_df[[variable_name1]] <- get(variable_name1)
        stat_per_pop_df[[variable_name2]] <- get(variable_name2)}

write.table(stat_per_pop_df, file = file.path("output_indels", "stat_per_pop.tsv"), sep = "\t")

long_het <- stat_per_pop_df %>% pivot_longer(-ID, names_to = c(".value", "Population"), names_pattern = "(Ho|He)_(.+)")
fis_df <- long_het %>% mutate(Fis = 1-(Ho/He))
shapiro.test(fis_df$Fis)
# 
# Shapiro-Wilk normality test
# 
# data:  fis_df$Fis
# W = 0.99068, p-value = 0.04646

fis_results <- fis_df %>% group_by(Population) %>% summarise(mean_Fis = mean(Fis, na.rm = TRUE),
               p_value = wilcox.test(na.omit(Fis), mu = 0, exact = FALSE)$p.value, W_stat = wilcox.test(na.omit(Fis), mu = 0)$statistic) %>%
               mutate(direction = ifelse(mean_Fis > 0, "Inbreeding", "Excess Heterozygosity"), sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value <0.05 ~ "*"), 
                      label_pos = ifelse(mean_Fis >= 0, mean_Fis + 0.003, mean_Fis - 0.003)) %>%
               arrange(p_value)
               
pdf(file.path("output_indels", "figures", "Fis_of_contributing_indels.pdf"))
ggplot(fis_results, aes(x = reorder(Population, mean_Fis), y = mean_Fis, fill = direction)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(y = label_pos, label = sig), hjust = 0.5, size  = 3.5, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Inbreeding"  = "salmon", "Excess Heterozygosity" = "turquoise3")) +
  labs(title = "Mean Inbreeding Coefficient (Fis) per Population", subtitle = "t-test | * p<0.05  ** p<0.01  *** p<0.001", x = "Population", y = "Mean Fis") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9), plot.title     = element_text(face = "bold"), plot.subtitle  = element_text(color = "grey40", size = 8), legend.position = "none")
dev.off()

long_df <- stat_per_pop_df %>% pivot_longer(-ID, names_to = c(".value", "Population"), names_pattern = "(Ho|He)_(.+)")

superpop <- data.frame( Population = c("ACB","ASW","ESN","GWD","LWK","MSL","YRI", "CEU","FIN","GBR","IBS","TSI", "CDX","CHB","CHS","JPT","KHV", "BEB","GIH","ITU","PJL","STU","CLM","MXL","PEL","PUR"), Superpop = c(rep("AFR", 7), rep("EUR", 5), rep("EAS", 5), rep("SAS", 5), rep("AMR", 4)))

pop_order <- superpop$Population

pdf(file.path("output_indels", "figures", "Ho-He_of contributing_indels.pdf"))
long_df %>% mutate(diff = Ho - He) %>% left_join(superpop, by = "Population") %>%
  mutate(Population = factor(Population, levels = pop_order)) %>%
  ggplot(aes(x = Population, y = diff, color = ID)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 0.8, linetype  = "solid") +
  annotate("rect", xmin = 0.5, xmax = 7.5, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = 0.07) +
  annotate("rect", xmin = 7.5, xmax = 12.5, ymin = -Inf, ymax = Inf, fill = "#0072B2", alpha = 0.07) +
  annotate("rect", xmin = 12.5, xmax = 17.5, ymin = -Inf, ymax = Inf, fill = "#009E73", alpha = 0.07) +
  annotate("rect", xmin = 17.5, xmax = 22.5, ymin = -Inf, ymax = Inf, fill = "#CC79A7", alpha = 0.07) +
  annotate("rect", xmin = 22.5, xmax = 26.5, ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = 0.07) +

  annotate("text", x = 4, y = Inf, label = "AFR", vjust = 2, fontface = "bold", size = 3.5, color = "#E69F00") +
  annotate("text", x = 10, y = Inf, label = "EUR", vjust = 2, fontface = "bold", size = 3.5, color = "#0072B2") +
  annotate("text", x = 15,   y = Inf, label = "EAS", vjust = 2, fontface = "bold", size = 3.5, color = "#009E73") +
  annotate("text", x = 20,   y = Inf, label = "SAS", vjust = 2, fontface = "bold", size = 3.5, color = "#CC79A7") +
  annotate("text", x = 24.5, y = Inf, label = "AMR", vjust = 2, fontface = "bold", size = 3.5, color = "#D55E00") +
  labs(title = "Ho - He per INDEL per Population", x = "Population", y = "Ho - He", color = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))
dev.off()

# Run pairwise Fst
pairwise_fst <- pairwise.neifst(genind2hierfstat(contrib_obj))

fst_matrix <- as.matrix(pairwise_fst)
superpop_map <- c(ACB="AFR", ASW="AFR", ESN="AFR", GWD="AFR", LWK="AFR", MSL="AFR", YRI="AFR",
  CEU="EUR", FIN="EUR", GBR="EUR", IBS="EUR", TSI="EUR", CDX="EAS", CHB="EAS", CHS="EAS", JPT="EAS",
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

pdf(file.path("output_indels", "figures", "Pairwise_Fst_of_contributing_indels.pdf"))
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

pdf(file.path("output_snps", "figures", "Fst_Ho_correlation_indels.pdf"))
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
