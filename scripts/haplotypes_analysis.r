#Set coordinates of interest
# region 1 => 7:114362480-114650436 / region 2 => 7:147017367-147323381 / region 3 => 7:148057724-148333434
# +-20kb / region 1 => 7:114342480-114670436 / region 2 => 7:146997367-147343381 / region 3 => 7:148037724-148353434

# Download raw file from repository
dir.create(file.path("data", "raw"), recursive = TRUE)
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/contr_phased.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "contr_phased.vcf.gz"))

output_dir <- file.path("output_haplo", "figures")
dir.create(output_dir, recursive = TRUE)

library(vcfR); library(ggplot2); library(reshape2); library(tidyr); library(tidyverse)
library(randomcoloR); library(adegenet); library(ade4); library(hierfstat); library(dplyr); library(stringdist); library(ape); library(plotly); library(GWLD)

vcf <- read.vcfR(file.path("data", "raw", "contr_phased.vcf.gz"))

#check for unique IDs
with(list(ids = getID(vcf)), length(unique(ids, incomparables = NA)) == length(ids)) 

genotypes <- extract.gt(vcf, element = "GT")
fix <- getFIX(vcf)
haplo <- matrix(NA, nrow = nrow(genotypes), ncol = ncol(genotypes)*2)
haplo1 <- apply(genotypes, 2, function(x) as.numeric(substr(x, 1, 1)))
haplo2 <- apply(genotypes, 2, function(x) as.numeric(substr(x, 3, 3)))
haplo[, seq(1, ncol(haplo), 2)] <- haplo1
haplo[, seq(2, ncol(haplo), 2)] <- haplo2

rownames(haplo) <- rownames(genotypes)
colnames(haplo) <- paste0(rep(colnames(genotypes), each = 2), c("_1", "_2"))

all_names <- data.frame(var = vcf@fix[, "ID"], pos = vcf@fix[, "POS"])

#subset for 3 different regions
region1_names <- all_names[all_names[,2] >= 114342480 & all_names[,2] <= 114670436, ]
region2_names <- all_names[all_names[, 2] >= 146997367 & all_names[, 2] <= 147343381, ]
region3_names <- all_names[all_names[, 2] >= 148037724 & all_names[, 2] <= 148353434, ]

regions_names <- list(region1 = region1_names, region2 = region2_names, region3 = region3_names)
#insert only IDs for contributing variants
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/contrib_ids.txt"
download.file(url, destfile = file.path("data", "raw", "contrib_ids.txt"))
contrib_ids <- readLines(file.path("data", "raw", "contrib_ids.txt"))

for(i in seq_along(regions_names))
    {region <- names(regions_names)[i]
     ids <- regions_names[[i]][,1]
     contr_haplo <- haplo[rownames(haplo) %in% ids, ]
     assign(paste0("contr_haplo_", region), contr_haplo)
     }

haplo_all_regions <- list( region1 = contr_haplo_region1, region2 = contr_haplo_region2, region3 = contr_haplo_region3)

# Population specific
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/pop.csv"
download.file(url, destfile = file.path("data", "raw", "pop.csv"))
pop <- read.table(file.path("data", "raw", "pop.csv"), header = FALSE, sep = ",")
pop <- pop[rep(1:nrow(pop), each = 2), ]
suf <- rep(c("_1", "_2"), times = nrow(pop)/2)
pop[, 2] <- paste0(pop[, 2], suf)
pop_unique <- unique(pop[,3])
for(i in seq_along(pop_unique))
    {pop_s <- pop_unique[[i]]
     samples <- pop[pop[, 3] == pop_unique[[i]], ]
     ids <- samples[, 2]
     for(y in seq_along(haplo_all_regions))
            {b <- haplo_all_regions[[y]]
             reg <- names(haplo_all_regions)[y]
             haplo_p <- b[, colnames(b) %in% ids]
             assign(paste0("contr_haplo_", reg, "_", pop_s), haplo_p)
            }
    }

reg_names <- ls(pattern = "^contr_haplo_region")
contr_reg <- mget(reg_names)
names(contr_reg) <- gsub("contr_haplo_", "", names(contr_reg))

for(i in seq_along(contr_reg))
    {hap <- contr_reg[[i]]
     n <- names(contr_reg)[i]
     string <- apply(hap, 2, paste0, collapse = "")
     n_string <- length(string)
     string <- table(string)/n_string
     string <- data.frame(haplo = names(string), freq = as.numeric(string))
     assign(paste0("string_", n), string)
     }

#create total matrices for 3 regions
strs_names <- ls(pattern = "^string_")
strs <- mget(strs_names)
names(strs) <- gsub("string_", "", names(strs))
strs_r <- list(region1 = string_region1, region2 = string_region2, region3 = string_region3)

for(i in seq_along(strs_r))
    {x <- names(strs_r)[i]
     df <- data.frame()
     z <- strs_r[[i]]
     if(grepl("region1*", x))
        {haplos <- unique(strs_r[[i]][, 1])
         df <- data.frame(haplotype = haplos)
         assign(paste0(x, "_df"), df)}
     else if (grepl("region2*", x))
        {haplos <- unique(strs_r[[i]][, 1])
         df <- data.frame(haplotype = haplos)
         assign(paste0(x, "_df"), df)}
     else if (grepl("region3*", x))
        {haplos <- unique(strs_r[[i]][, 1])
         df <- data.frame(haplotype = haplos)
         assign(paste0(x, "_df"), df)
        }
    }


for(i in seq_along(strs))
    {v <- strs[[i]]
     d <- names(strs)[i]
     if(grepl("region1", d))
        {region1_df[[d]] <- v[,2][ match(region1_df[,1], v[, 1])]}
    else if (grepl("region2", d))
        {region2_df[[d]] <- v[,2][ match(region2_df[,1], v[, 1])]}
    else if (grepl("region3", d))
        {region3_df[[d]] <- v[,2][ match(region3_df[,1], v[, 1])]}
    }



superpop_map <- c(
  ACB="AFR", ASW="AFR", ESN="AFR", GWD="AFR", LWK="AFR", MSL="AFR", YRI="AFR",
  CEU="EUR", FIN="EUR", GBR="EUR", IBS="EUR", TSI="EUR",
  CDX="EAS", CHB="EAS", CHS="EAS", JPT="EAS", KHV="EAS",
  CLM="AMR", MXL="AMR", PEL="AMR", PUR="AMR",
  BEB="SAS", GIH="SAS", ITU="SAS", PJL="SAS", STU="SAS")


dfs <- list(region1 = region1_df, region2 = region2_df, region3 = region3_df)
hap_order_per_region <- list()
for(i in seq_along(dfs))
{region_df <- dfs[[i]]
reg <- names(dfs)[i]

plot_df <- region_df %>%
  select(haplotype, matches("_")) %>%
  pivot_longer(cols = where(is.numeric), names_to = "population", values_to = "freq") %>%
  mutate(population = gsub("^.*_", "", population), superpop   = superpop_map[population])

global_freq <- reg

#freq
global_hap_order <- region_df %>% arrange(desc(.data[[global_freq]])) %>% pull(haplotype)
hap_labels <- setNames(as.character(seq_along(global_hap_order)), global_hap_order)
hap_order_per_region[[i]] <- global_hap_order
plot_df <- plot_df %>% mutate(hap_label = hap_labels[as.character(haplotype)])

all_hap_labels <- as.character(seq_along(global_hap_order))
n_haps_global  <- length(all_hap_labels)

set.seed(42)
global_cols <- distinctColorPalette(n_haps_global)
names(global_cols) <- all_hap_labels

plot_df <- plot_df %>%
  mutate(superpop = factor(superpop, levels = c("AFR", "EUR", "EAS", "AMR", "SAS")))

superpops <- c("AFR", "EUR", "EAS", "AMR", "SAS")
plot_list <- list()

for (sp in superpops) 
  {sp_df <- plot_df %>% filter(superpop == sp)
  
  hap_order_sp <- sp_df %>% group_by(hap_label) %>% summarise(mean_freq = mean(freq, na.rm = TRUE)) %>%
    arrange(mean_freq) %>% pull(hap_label)
  
  sp_df <- sp_df %>% mutate(hap_label = factor(hap_label, levels = hap_order_sp))

  sp_df <- sp_df %>% arrange(population, desc(hap_label)) %>% group_by(population) %>%
    mutate(cum_freq  = cumsum(freq), label_pos = cum_freq - freq / 2, label     = ifelse(freq >= 0.05, as.character(hap_label), NA)) %>% ungroup()
  
  p <- ggplot(sp_df, aes(x = population, y = freq, fill = hap_label)) +
    geom_bar(stat = "identity", width = 0.85, colour = "white", linewidth = 0.2) +
    geom_text(aes(y = label_pos, label = label),
      size = 3, fontface = "bold", colour = "grey15", na.rm = TRUE) +
    scale_fill_manual(values = global_cols, name = "Haplotype") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits = c(0, 1.01)) +
    labs(x = "Population", y = "Haplotype Frequency", title = paste0(sp)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), legend.position = "none", plot.title = element_text(face = "bold"))
  
  plot_list[[sp]] <- p
  
  ggsave(filename = file.path("output_haplo", "figures", paste0(reg, sp, "_haplotypes.png")), plot = p, width = 6, height = 5)
 }}


pop_clean <- read.table(file.path("data", "raw", "pop.csv"), header = FALSE, sep = ",")
total_df <- data.frame(row.names = pop_clean[,2])

pdf(file.path("output_haplo", "figures", "haplo_analysis.pdf"))
#haplo_all_regions list
for(i in seq_along(haplo_all_regions))
  {contr_haplo <- haplo_all_regions[[i]]
   reg <- names(haplo_all_regions)[i]
   hap <- t(contr_haplo)
   hap1 <- hap[grepl("_1$", rownames(hap)), ]
   rownames(hap1) <- gsub("_1$", "", rownames(hap1))
   hap2 <- hap[grepl("_2$", rownames(hap)), ]
   rownames(hap2) <- gsub("_2$", "", rownames(hap2))
   hap1_string <- apply(hap1, 1, paste, collapse = "")
   hap2_string <- apply(hap2, 1, paste, collapse = "")
   genot_df <- data.frame(hap = paste(hap1_string, hap2_string, sep = "/"), row.names = rownames(hap1))
   total_df <- cbind(total_df, genot_df)
   genind <- df2genind(X = genot_df, pop = pop_clean[, 3], sep = "/")
   assign(paste0("genind_", reg), genind)
   genpop <- genind2genpop(genind)
   assign(paste0("genpop_", reg), genpop)}

colnames(total_df) <- c("hap1", "hap2", "hap3")
genind_total <- df2genind(X = total_df, pop = pop_clean[, 3], sep = "/")


#statistics
geninds <- list(region1 = genind_region1, region2 = genind_region2, region3 = genind_region3)
stat_perloc_df <- data.frame()
stat_overall_df <- data.frame()
heteroz_df <- data.frame()
for(i in seq_along(geninds))
  {genind_i <- geninds[[i]]
   genind_i@all.names <- lapply(genind_i@all.names, function(x) as.character(1:length(x)))
   hierf <- genind2hierfstat(genind_i)
   statistics <- basic.stats(hierf)
   temp1 <- as.data.frame(statistics$perloc["hap", ])
   stat_perloc_df <- rbind(stat_perloc_df, temp1)
   temp2 <- as.data.frame(t(statistics$overall))
   stat_overall_df <- rbind(stat_overall_df, temp2)
   temp3 <- as.data.frame(t(statistics$Ho["hap", ]))
   temp4 <- as.data.frame(t(statistics$Hs["hap", ]))
   heteroz_df <- rbind(heteroz_df, temp3)
   heteroz_df <- rbind(heteroz_df, temp4)
}

rownames(stat_overall_df) <- names(geninds)
rownames(stat_perloc_df) <- names(geninds)
row.names(heteroz_df) <- c("Ho_region1", "Hs_region1", "Ho_region2", "Hs_region2", "Ho_region3", "Hs_region3")

#pairwise Fst
pairwise_fst <- function(genind)
{pop <- levels(pop(genind)); n <- length(pop)
 fst_matrix <- matrix(NA, n, n, dimnames = list(pop, pop))
 for(i in 1:(n-1))
 {for(j in (i+1):n)
    {pair <- genind[pop(genind) %in% c(pop[i], pop[j]), ]
     pop(pair) <- droplevels(pop(pair))
     pair@all.names <- lapply (pair@all.names, function(x) as.character(seq_along(x)))
     hierf <- genind2hierfstat(pair)
     fst_matrix2 <- pairwise.WCfst(hierf)
     fst_matrix[i, j] <- fst_matrix2[1, 2]; fst_matrix[j, i] <- fst_matrix2[1, 2]
 }}
 return(fst_matrix)
 }

pairwise_reg1 <- pairwise_fst(genind_region1)
pairwise_reg2 <- pairwise_fst(genind_region2)
pairwise_reg3 <- pairwise_fst(genind_region3)

plot_fst_heatmap <- function(fst_matrix, title, superpop_map)
{pop <- rownames(fst_matrix); superpop <- superpop_map[pop]
 order_i <- pop[order(superpop)]
 fst_matrix <- fst_matrix[order_i, order_i]
 fst_long <- reshape2::melt(fst_matrix, na.rm = FALSE); colnames(fst_long) <- c("Pop1", "Pop2", "Fst")
 fst_long$Pop1 <- factor(fst_long$Pop1, levels = order_i)
 fst_long$Pop2 <- factor(fst_long$Pop2, levels = order_i)
 superpop_df <- data.frame(pop = factor(order_i, levels = order_i), superpop = superpop_map[order_i])
 colors_superp <- c("AFR" = "#E69F00", "AMR" = "#56B4E9", "EUR" = "#F0E442", "EAS" = "#009E73", "SAS" = "#CC79A7")
 heatmap <- ggplot(fst_long, aes(x = Pop1, y = Pop2, fill = Fst)) +
                    geom_tile(color = "white") +
                    scale_fill_gradient(low = "lightpink", high = "darkred", na.value = "grey90", name = expression(F[ST])) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                          axis.text.y = element_text(size = 7), panel.grid = element_blank(), 
                          plot.title = element_text(hjust = 0.05), legend.position = "right")
print(heatmap)
}

p1 <- plot_fst_heatmap(pairwise_reg1, "Locus 1", superpop_map)
p2 <- plot_fst_heatmap(pairwise_reg2, "Locus 2", superpop_map)
p3 <- plot_fst_heatmap(pairwise_reg3, "Locus 3", superpop_map)

# Fst to Ho visualization + Spearman correlation
mean_fst <- function(fst_matrix, region_name)
  {diag(fst_matrix )<- NA
   data.frame(Population = rownames(fst_matrix), mean_Fst = rowMeans(fst_matrix, na.rm = TRUE), Region = region_name)
}

fst_long <- bind_rows(mean_fst(pairwise_reg1, "region1"), mean_fst(pairwise_reg2, "region2"), mean_fst(pairwise_reg3, "region3"))

heteroz_long <- heteroz_df |> tibble::rownames_to_column("stat_region") |>
  pivot_longer(-stat_region, names_to = "Population", values_to = "Value") |>
  separate(stat_region, into = c("Stat", "Region"), sep = "_", extra = "merge") |>
  mutate(Superpopulation = superpop_map[Population])

ho_long <- heteroz_long |> filter(Stat == "Ho") |> select(Population, Region, Ho = Value, Superpopulation)
correlation <- left_join(ho_long, fst_long, by = c("Population", "Region"))
correl_results <- correlation |> group_by(Region) |>
    summarise(rho = cor(Ho, mean_Fst, method = "spearman", use = "complete.obs"), 
    p_val = cor.test(Ho, mean_Fst, method = "spearman", exact = FALSE)$p.value, .groups = "drop") |> 
    mutate(significant = ifelse(p_val <- 0.05, "*", "ns"))

ggplot(correlation, aes(x = Ho, y = mean_Fst, color = Superpopulation, label = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed", linewidth = 0.7) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  scale_color_manual(values = c("AFR" = "#E69D00", "AMR" = "#56B4E9", "EUR" = "#F0E442", "EAS" = "#009E73", "SAS" = "#CC79A7")) + 
  facet_wrap(~ Region, scale = "free") +
  labs(title = expression("Correlation Between Ho and mean Fst"),
       x = expression("Observed Heterozygosity (Ho)"), y = expression("Mean pairwise Fst")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position  = "bottom", strip.text = element_text(face = "bold"))


#PCA

Col1 <- rainbow(length(levels(pop(genind))))
myPal <- colorRampPalette(c("blue","gold","red"))
#total analysis - UNFILTERED
X <- tab(genind_total, freq=TRUE, NA.method = "mean")
pca1 <- dudi.pca(genind_total, scannf = FALSE, nf = 3)
s.class(pca1$li, fac = pop(genind_total), xax = 1, yax = 2, col = transp(Col1, 0.6), axesell = FALSE, cstar = 0, cpoint = 3, grid = FALSE)

#Filter outliers
apply(pca1$li, 2, range)
#       Axis1     Axis2      Axis3
#[1,] -3.196912 -3.989578 -28.919199
#[2,] 60.921483 79.421274   8.557886
rownames(pca1$li)[which(pca1$li[,1] > 20)]
#[1] "HG02944"
rownames(pca1$li)[which(pca1$li[,2] > 20)]
#[1] "HG02944"
rownames(pca1$li)[which(pca1$li[,3] < -10)]
#[1] "NA19144"
outliers <- c("HG02944", "NA19144")
genind_clean <- genind_total[!indNames(genind_total) %in% outliers, ]
pca2 <- dudi.pca(genind_clean, scannf = FALSE, nf = 3)
s.class(pca2$li, fac = pop(genind_clean), xax = 1, yax = 2, col = transp(Col1, 0.6), axesell = FALSE, cstar = 0, cpoint = 3, grid = FALSE)

#Second round of filtering
apply(pca2$li, 2, range)
#Axis1     Axis2     Axis3
#[1,] -3.936060 -8.297135 -4.763799
#[2,]  6.139918 25.062148 75.828089
rownames(pca2$li)[which(pca2$li[,2] > 20)]
#[1] "HG02282"
rownames(pca2$li)[which(pca2$li[,3] > 30)]
#[1] "HG02282"
outliers2 <- c("HG02282")
genind_clean2 <- genind_clean[!indNames(genind_clean) %in% outliers2, ]
pca3 <- dudi.pca(genind_clean2, scannf = FALSE, nf = 3)
s.class(pca3$li, pop(genind_clean2),xax=1,yax=2, col=Col1)
title("PCA of diversity \ axes 1-2")
s.class(pca3$li, fac = pop(genind_clean2), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)
colorplot(pca3$li[c(1,2)], pca3$li, transp=TRUE, cex=1.5, xlab="PC 1", ylab="PC 2")
title("PCA of human genes \ axes 1-2")
abline(v=0,h=0,col="grey", lty=2)

#3 individuals removed in total

genpop_clean2 <- genind2genpop(genind_clean2)
ca <- dudi.coa(tab(genpop_clean2), scannf = FALSE, nf = 2)

barplot(ca$eig, main = "CA eigenvalues ", col = heat.colors(length(ca$eig)))
s.label(ca$li, sub = "CA 1,2", csub = 2)
add.scatter.eig(ca$eig, nf = 3, xax = 1, yax = 3, posi = "bottomright")

grp <- find.clusters(genind_clean2, max.n.clust = 100, n.pca = 150, n.clust = 5)
#150 PCs
#5 clusters
dapc1 <- dapc(genind_clean2, grp$grp, , scannf = FALSE, n.pca = 100, n.da = 4)
#100 PCs - 4 discr

myCol <- c("orange", "red", "blue", "green", "pink")
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  par(mgp = c(1.2, 0.5, 0))
  plot(temp, col=rep(c("black", "lightgray"),
       c(dapc1$n.pca,1000)), ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance",
       cex=0.5, pch=10, type="h", lwd=2, cex.lab = 0.7, cex.axis = 0.7) }

scatter(dapc1,xax=1,yax=3, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0, col=Col1,solid=.4, cex=1.5, clab=0, mstree=TRUE, scree.da=FALSE
        , posi.pca="topright",leg=TRUE,txt.leg=paste("Cluster",1:5))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,3], pch=4, cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,3], pch=4, cex=3, lwd=2, col=myCol)
add.scatter(myInset(), posi="topleft", inset=c(-0.05,-0.1), ratio=.10, bg=transp("white"))
scatter(dapc1,3,3, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

dapc2 <- dapc(genind_clean2, n.da = 4, n.pca = 150)
#150 PCs

optim_pcs <- optim.a.score(dapc2)

dapc3 <- dapc(genind_clean2, n.pca = 36, n.da = 4)
scatter(dapc3,xax=1,yax=2, col=transp(myPal(8)), scree.da=FALSE, cell=1.5, cex=1, bg="white",cstar=0, legend=TRUE)
compoplot(dapc3, cleg = .5, posi = list(x = 0, y = 1.2, lab = pop(genind_clean2)))

contrib <- loadingplot(dapc3$var.contr, axis=2, thres=0.04, lab.jitter=1)
contr_haplo <- contrib$var.names
loadings <- dapc3$var.contr[, 2]
loadingplot(dapc3$var.contr, axis = 2, thres = 0.04, lab.jitter = 1, lab = " ")

df2 <- data.frame(haplotype = names(loadings), loading = loadings, index = seq_along(loadings))

haplo_palette <- c(
  "hap1.0011011"= "seagreen", 
  "hap1.0111011" = "palegreen3",
  "hap1.0000000" = "aquamarine",
  "hap1.0100000" = "darkolivegreen1",
  "hap2.00000000010000000000000000000000100" = "magenta4",
  "hap2.00000000000000000000000000000000000" = "magenta2",
  "hap2.10101111101111111111111111111111011" = "purple3",
  "hap2.10000100010000000000000000000000100" = "lightpink3",
  "hap3.0000000000000000000000000000000001" = "royalblue4", 
  "hap3.0000000000000000000000000000000010" = "skyblue", 
  "hap3.0000000000000000000000000000000000" = "royalblue"
)

freq_contr_haplo <- `[`(makefreq(genpop_clean2, quiet = TRUE), , contr_haplo)
freq_contr_long <- melt(freq_contr_haplo, varnames = c("pop", "hap"), value.name = "freq")


ggplot(freq_contr_long, aes(x = pop, y = freq, fill = hap)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ gsub("\\..*", "", hap), scales = "free_x") +
  scale_fill_manual(values = haplo_palette) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold"),) +
    labs(title = "Contributing Haplotype Frequencies per Population", x = "Population", y = "Frequency", fill  = "Haplotype")

freq_contr_long_loc1 <- subset(freq_contr_long, grepl("hap1", hap))
pop_order <- names(superpop_map)
freq_contr_long_loc1 <- freq_contr_long_loc1 %>% mutate(superpop = superpop_map[pop], pop = factor(pop, levels = pop_order))

ggplot(freq_contr_long_loc1, aes(x = pop, y = freq, group = hap, color = hap)) +
   geom_line(size = 1) +
   geom_point(size = 2) +
   facet_wrap(~hap, scales = "free_y") +
   scale_color_manual(values = haplo_palette) +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
   	   annotate("rect", xmin = 0.5,  xmax = 7.5, ymin = -Inf, ymax = Inf, fill = "#E69F00", alpha = 0.07) +
       annotate("rect", xmin = 7.5, xmax = 12.5, ymin = -Inf, ymax = Inf, fill = "#0072B2", alpha = 0.07) +
       annotate("rect", xmin = 12.5, xmax = 17.5, ymin = -Inf, ymax = Inf, fill = "#009E73", alpha = 0.07) +
       annotate("rect", xmin = 17.5, xmax = 21.5, ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = 0.07) +
       annotate("rect", xmin = 21.5, xmax = 26.5, ymin = -Inf, ymax = Inf, fill = "#CC79A7", alpha = 0.07) +
      annotate("text", x = 4, y = Inf, label = "AFR", vjust = 2, fontface = "bold", size = 3.5, color = "#E69F00") +
      annotate("text", x = 10, y = Inf, label = "EUR", vjust = 2, fontface = "bold", size = 3.5, color = "#0072B2") +
      annotate("text", x = 15, y = Inf, label = "EAS", vjust = 2, fontface = "bold", size = 3.5, color = "#009E73") +
      annotate("text", x = 20, y = Inf, label = "AMR", vjust = 2, fontface = "bold", size = 3.5, color = "#D55E00") +
      annotate("text", x = 24, y = Inf, label = "SAS", vjust = 2, fontface = "bold", size = 3.5, color = "#CC79A7") +
  
      facet_wrap(~hap, scales = "free_y") +
      labs(title = "Contributing Haplotype Frequencies - Locus 1", subtitle = "Organised by superpopulation",
           x = "Population", y = "Frequency", fill = "Superpopulation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(color = "grey40"), legend.position = "bottom")


alleles_df <- data.frame()
vcf <- as.data.frame(getFIX(vcf))
al <- data.frame(VAR = vcf$ID, POS = vcf$POS, REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)
alleles_df <- rbind(alleles_df, al)
alleles_df <- alleles_df[order(alleles_df[, 2]), ]

al_con <- alleles_df[alleles_df$VAR %in% contrib_ids, ]
r1_df <- al_con[al_con[, 2] >= 114342480 & al_con[,2] <= 114670436, ]
r1_df <- separate(r1_df, ALT, into = c("ALT", "ALT2"), sep = ",", fill = "right")
r2_df <- al_con[al_con[, 2] >= 146997367 & al_con[, 2] <= 147343381, ]
r2_df <- separate(r2_df, ALT, into = c("ALT", "ALT2"), sep = ",", fill = "right")
r3_df <- al_con[al_con[, 2] >= 148037724 & al_con[, 2] <= 148353434, ]
r3_df <- separate(r3_df, ALT, into = c("ALT", "ALT2"), sep = ",", fill = "right")
r1_df$ALT <- c("C", "G----", "G-", "C", "G-", "TA", "C")
r2_df$REF[1] <- "A--"; r2_df$ALT[1] <- "AT-"
r2_df$ALT2[is.na(r2_df$ALT2)] <- "-"
r3_df$ALT[11] <- "T--"; r3_df$ALT[14] <- "T---"
r3_df$ALT2[is.na(r3_df$ALT2)] <- "-"; r3_df$ALT[29] <- "A-"; r3_df$ALT[30] <- "C--"; r3_df$ALT[31] <- "A--"
r3_df$REF[32] <- "A-"; r3_df$ALT[33] <- "C--"; r3_df$ALT[34] <- "G-"

r_dfs <- list(hap1 = r1_df, hap2 = r2_df, hap3 = r3_df)
contr_haplo_list <- split(gsub(".*\\.", "", contr_haplo), gsub("\\..*", "", contr_haplo))
total_haplo_seq <- list()

for(r in seq_along(r_dfs))
  {r_df <- r_dfs[[r]]
   haplos_i <- contr_haplo_list[[r]]
   total_haplo_i <- c()
   for(l in seq_along(haplos_i))
      {haplo0 <- strsplit((haplos_i), "")[[l]]
       mat <- matrix(ncol = length(haplo0), nrow = 2)
       mat[1, ] <- r_df[,1]
       mat[2, ] <- haplo0
       for(y in 1:ncol(mat))
          {id <- mat[1, y]
           ref <- r_df$REF[r_df$VAR == id][1]
           alt <- r_df$ALT[r_df$VAR == id][1]
           alt2 <- r_df$ALT2[r_df$VAR == id][1]
           if(mat[2, y] == "0")
              {mat[2, y] <- ref}
           else if(mat[2, y] == "1")
              {mat[2, y] <- alt}
           else if(mat[2, y] == "2")
              {mat[2, y] <- alt2}
          }
    temp_matrix <- mat[2, ]
    temp_string <- paste(temp_matrix, collapse = "")
    total_haplo_i <- c(total_haplo_i, temp_string)
    }
  total_haplo_seq[[r]] <- total_haplo_i
} 

haplo_seq_df <- data.frame(region = rep(seq_along(total_haplo_seq), lengths(total_haplo_seq)), haplotype = unlist(total_haplo_seq), row.names = NULL)
#edit for region1 (coexisting variants rs71999274 & rs7799652)

haplo_seq_df$haplotype[1] <- "TGTTTAG-CTAC"
haplo_seq_df$haplotype[2] <- "TG----G-CTAC"
haplo_seq_df$haplotype[3] <- "TGTTTAGATGAA"
haplo_seq_df$haplotype[4] <- "TG----GATGAA"

haplo_long <- haplo_long <- haplo_seq_df %>% rowwise() %>% mutate(pos  = list(seq_len(nchar(haplotype))),
         char = list(strsplit(haplotype, "")[[1]])) %>% unnest(c(pos, char)) %>% ungroup()

haplo_palette2 <- setNames(c("seagreen", "palegreen3", "aquamarine", "darkolivegreen1", "magenta4",
                   "magenta2", "purple3", "lightpink3", "royalblue4", "skyblue", "royalblue"), haplo_seq_df$haplotype)

ggplot(haplo_long, aes(x = pos, y = haplotype, fill = haplotype)) +
  geom_tile(colour = "grey8", linewidth = 0) +                    
  geom_text(aes(label = char), colour = "grey8", family = "Courier", fontface = "bold", size = 3) +
  scale_fill_manual(values = haplo_palette2) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    strip.text = element_blank(), legend.position = "none")

# PCoA, Hamming and visualization for all haplotypes
all_haplo_list <- genind_clean2@all.names
all_haplo_list <- lapply(seq_along(all_haplo_list), function(r)
  {all_haplo_list[[r]][match(hap_order_per_region[[r]], all_haplo_list[[r]])]})
all_haplo_seq <- list()

for(r in seq_along(r_dfs))
  {r_df <- r_dfs[[r]]
   haplos_i <- all_haplo_list[[r]]
   total_haplo_i <- c()
   for(l in seq_along(haplos_i))
      {haplo0 <- strsplit((haplos_i), "")[[l]]
       mat <- matrix(ncol = length(haplo0), nrow = 2)
       mat[1, ] <- r_df[,1]
       mat[2, ] <- haplo0
       for(y in 1:ncol(mat))
        {id <- mat[1, y]
         ref <- r_df$REF[r_df$VAR == id][1]
         alt <- r_df$ALT[r_df$VAR == id][1]
         alt2 <- r_df$ALT2[r_df$VAR == id][1]
         if(mat[2, y] == "0")
            {mat[2, y] <- ref}
         else if(mat[2, y] == "1")
            {mat[2, y] <- alt}
         else if(mat[2, y] == "2")
            {mat[2, y] <- alt2}
    }
    temp_matrix <- mat[2, ]
    temp_string <- paste(temp_matrix, collapse = "")
    total_haplo_i <- c(total_haplo_i, temp_string)
}
  all_haplo_seq[[r]] <- total_haplo_i
} 

all_haplo_seq_df <- data.frame(region = rep(seq_along(all_haplo_seq), lengths(all_haplo_seq)), haplotype = unlist(all_haplo_seq), row.names = NULL)

# Hamming Distances

hamming2 <- lapply(seq_along(all_haplo_seq), function(i)
  {seq <- all_haplo_seq[[i]]
  mat_i <- stringdistmatrix(seq, seq, method = "hamming")
  rownames(mat_i) <- colnames(mat_i) <- paste0("hap", seq_along(seq))
  mat_i})


pcoa_per_region <- lapply(hamming2, function(mat) {
  pcoa_res <- ape::pcoa(mat)
  pcoa_res
})

coords_per_region <- lapply(pcoa_per_region, function(p) p$vectors[, 1:3])

var_per_region <- lapply(pcoa_per_region, function(p) {
  eigenvalues <- p$values$Eigenvalues
  eigenvalues <- pmax(eigenvalues, 0)
  (eigenvalues / sum(eigenvalues)) * 100
})

contributing_per_region <- lapply(1:3, function(i) 
 {contrib_seqs <- haplo_seq_df$haplotype[haplo_seq_df$region == i]
  contrib_idx <- which(all_haplo_seq[[i]] %in% contrib_seqs)
  paste0("hap", contrib_idx)
 })

colors_l1.2 <- setNames(c( "hap1" = "seagreen", "hap6" = "palegreen3", "hap9" = "aquamarine", "hap23" = "darkolivegreen1"), contributing_per_region[[1]])
colors_l2.2 <- setNames(c("hap7" = "magenta4", "hap9" = "magenta2", "hap81" = "purple3", "hap109" = "lightpink3"), contributing_per_region[[2]])
colors_l3.2 <- setNames(c("hap1" = "royalblue4", "hap2" = "skyblue", "hap3" = "royalblue"), contributing_per_region[[3]])

color_groups <- list(colors_l1.2, colors_l2.2, colors_l3.2)

# Contributing haplotypes per region from haplo_seq_df
for (i in 1:3) {
  coords_i <- coords_per_region[[i]]
  var_i <- var_per_region[[i]]
  hap_names <- paste0("hap", seq_len(nrow(coords_i)))
  contributing <- contributing_per_region[[i]]
  contrib_colors <- color_groups[[i]]
  point_colors <- rep("gray40", length(hap_names))
  names(point_colors) <- hap_names
  point_colors[contributing] <- contrib_colors[contributing]
  point_colors_hex <- sapply(point_colors, col2hex)
  
  fig <- plot_ly()
  
  for (j in seq_len(nrow(coords_i))) {
    fig <- fig %>% add_trace(type = "scatter3d", mode = "markers", x = coords_i[j, 1], 
           y = coords_i[j, 2], z = coords_i[j, 3], text = hap_names[j], hoverinfo = "text",
           marker = list(color = point_colors_hex[j], 
           size = ifelse(hap_names[j] %in% contributing, 7, 4)), showlegend = FALSE)
  }

  contrib_idx <- which(hap_names %in% contributing)
  for (j in contrib_idx) {
    col_rgb <- col2rgb(point_colors[j])
    fill_color <- sprintf("rgba(%d,%d,%d,0.7)", col_rgb[1], col_rgb[2], col_rgb[3])
    fig <- fig %>% add_trace(type = "scatter3d", mode = "markers",
      x = coords_i[j, 1], y = coords_i[j, 2], z = coords_i[j, 3], text = hap_names[j], hoverinfo = "text", 
      marker = list(color = fill_color, size = 12, width = 0), showlegend = FALSE)
  } 

  fig <- fig %>% layout(
    title = paste("Locus", i),
    scene = list(
      xaxis = list(title = paste0("PCoA1 (", round(var_i[1], 1), "%)")),
      yaxis = list(title = paste0("PCoA2 (", round(var_i[2], 1), "%)")),
      zaxis = list(title = paste0("PCoA3 (", round(var_i[3], 1), "%)"))))
  print(fig)

}

for (i in 1:3) {
  coords_i <- coords_per_region[[i]]
  var_i <- var_per_region[[i]]
  hap_names <- paste0("hap", seq_len(nrow(coords_i)))
  contributing <- contributing_per_region[[i]]
  contrib_colors <- color_groups[[i]]
  
  point_colors <- rep("gray40", length(hap_names))
  names(point_colors) <- hap_names
  point_colors[contributing] <- contrib_colors[contributing]
  point_colors_hex <- sapply(point_colors, col2hex)
  
  plot(coords_i[, 1], coords_i[, 2],
       xlab = paste0("PCoA1 (", round(var_i[1], 1), "%)"),
       ylab = paste0("PCoA2 (", round(var_i[2], 1), "%)"),
       main = paste0("Locus ", i, " (", round(sum(var_i[1:2]), 1), "% variance)"),
       pch = 19, col = point_colors_hex, cex = 0.8,
       xlim = range(coords_i[, 1]) * 1.3,
       ylim = range(coords_i[, 2]) * 1.3)
  
  # Highlight contributing haplotypes
  contrib_idx <- which(hap_names %in% contributing)
  points(coords_i[contrib_idx, 1], coords_i[contrib_idx, 2],
         pch = 19, col = point_colors_hex[contrib_idx], cex = 1.5)
  
  # Add halo around contributing
  points(coords_i[contrib_idx, 1], coords_i[contrib_idx, 2],
         pch = 21, col = point_colors_hex[contrib_idx],
         bg = NA, cex = 2.5, lwd = 1.5)
  
  abline(h = 0, v = 0, lty = 2, col = "grey70")
}


#Entropy Linkage Disequilibrium - GWLD

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
geno_region1 <- geno_num[ ,colnames(geno_num) %in% region1_names[, 1]]
geno_region2 <- geno_num[ ,colnames(geno_num) %in% region2_names[, 1]]
geno_region3 <- geno_num[ ,colnames(geno_num) %in% region3_names[, 1]]

Info <- data.frame(CHROM = vcf[, 1], POS   = as.numeric(vcf[, 2]), ID    = vcf[, 3])
mycols <- colorRampPalette(c("#ffe6e6", "#ff0000", "#990000"))(100)

rmi_list <- list()


for(i in seq_along(pop_unique))
    {pop_i <- pop_unique[[i]]
     ids <- pop_clean[pop_clean[, 3] == pop_i, ]
     geno_num_i <- geno_num[row.names(geno_num) %in% ids[, 2], ]
     rmi <- GWLD(geno_num_i, method = "RMI", cores = 4)
     rmi_list[[pop_i]] <- GWLD(geno_num_i, method = "RMI", cores = 4)
     p <- HeatMap(geno_num_i, method = "RMI", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 4, color = mycols, label.size = 8)

grid::grid.text(pop_i, x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 14, fontface = "bold", col = "black"))
}

dev.off()
