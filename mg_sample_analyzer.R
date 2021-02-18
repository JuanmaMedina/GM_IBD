library("tidyverse")
library("dunn.test")

setwd("../Desktop/IBBTEC/articulos/paper_IBD/methods/MGX/")

########################################
#### GENE EXPRESSION BETWEEN GROUPS ####
########################################

## FUNCTIONS ##

gene_reader <- function(x){
  # This function reads DIAMOND merged output from FG_reader.R
  x <- read.table(x)

  colnames(x) <- c("sample","gene_ID","reads","patient","status")

  x$status[x$status=="nonIBD"] <- "H"                                           # Replace "nonIBD" label
  x <- subset(x, x[ , 3] > 10)                                                  # Remove samples with <10 reads

  return(x)

}


gene_distr <- function(x){
  # This function calculates distributions and statistics on the genes
  
  options(warn=-1)                                                              # turn off tibble outdated warning
  
  s <- as.data.frame(                                                           # statistics of each group
    group_by(x, status) %>%
      summarise(
        count = n(),
        mean = mean(reads, na.rm = T),
        sd = sd(reads, na.rm = T),
        median = median(reads, na.rm = T),
        IQR = IQR(reads, na.rm = T)
      ))
  
  print(s)
  
  options(warn=0)                                                               # turn on global warnings again
  
  wct <- pairwise.wilcox.test(x$reads, x$status, p.adjust.method = "BH")      # Pairwise Wilcoxon RS test (with BH)
  print(wct)
  
  dut <- dunn.test(x$reads, x$status, method="bh", altp = T)                    # Dunn test (with BH)
  print(dut)                                                                    # altp=T means pvals, not zscores
}


den_plot <- function(x,title,a=0,b=4){
  x <- ggplot(x, aes(x = reads * 10^-3)) +                                      # log2(reads+1)
    geom_density(aes(fill = status), alpha=0.6) +                               # control line width with "size="
    theme_minimal()+
    theme (axis.title.x=element_blank()) + 
    ggtitle(title) +
    ylab("Density") +
    coord_cartesian(xlim=c(a, b)) +
    xlab("Abundance")
  
  return(x)
}

# box_plot <- function(x, title){
#   x <- ggplot(x, aes(x = status, y = reads/1000, fill = status)) +
#     geom_boxplot(coef = 18) + 
#     theme_minimal() +
#     theme(axis.title.x=element_blank(),
#         axis.text.x = element_text(colour="black", size=12, face="italic")) +
#   # scale_y_continuous(breaks = seq(0, 20000, by = 2500)) +
#     coord_cartesian(ylim = c(0, 11.5)) +
#     ylab("Abundance") +
#     ggtitle(title)
#   
#   return(x)
# }

## TDCD (PROPIONATE KINASE) ANALYSIS ## 
tdcd <- gene_reader(x="tdcd_all.tsv")
gene_distr(tdcd)
den_plot(tdcd,expression(italic("tdcD")),a=2,b=15)
# box_plot(tdcd,expression(italic("tdcD")))

## PDUW (PROPIONATE KINASE) ANALYSIS ## 
pduw <- gene_reader(x="pduw_all.tsv")
gene_distr(pduw)
den_plot(pduw,expression(italic("pduW")),a=1,b=11.5)
# box_plot(pduw,expression(italic("tdcD")))

## PCT (PROPIONATE COA TRF) ANALYSIS ## 
pct <- gene_reader(x="pct_all.tsv")
gene_distr(pct)
den_plot(pct,expression(italic("pct")),b=1)

## PRPE (PROPIONATE COA LIGASE) ANALYSIS ## 
prpe <- gene_reader(x="prpe_all.tsv")
gene_distr(prpe)
den_plot(prpe,expression(italic("prpE")),b=1)

## ACKA (ACETATE KINASE) ANALYSIS ## 
acka <- gene_reader(x="acka_all.tsv")
gene_distr(acka)
den_plot(acka,expression(italic("ackA")),a=1.5,b=5)

## SUM (PDUW + TDCD)                                                            
kin2 <- gene_reader(x="propk2_all.tsv")
gene_distr(kin2)
den_plot(kin2,expression(italic("tdcD + pduW")),a=0.2,b=2.2)


# Relative abundance of the 4 genes
prop4 <- bind_rows(list(tdcD=tdcd, pduW=pduw, 
                          pct=pct, prpE=prpe), .id = "gene")
prop4$gene <- factor(prop4$gene , levels=c("tdcD", "pduW", "pct", "prpE"))      # reorder categories in plot with factor

# A: Boxplot variant
ggplot(prop4, aes(x = gene, y = reads/1000)) +
  geom_boxplot(coef = 25) +                                                     # remove dots and include outliers in whiskers 
  guides(fill = FALSE) +                                                        # https://stackoverflow.com/questions/34602872/
  theme_minimal() +                                                             # no criteria found, fine-tune manually
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(colour="black", size=12, face="italic")) +
  # scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  coord_cartesian(ylim = c(0, 11.5)) +
  ylab("Abundance")

# Split A by condition
ggplot(prop4, aes(x = gene, y = reads/1000, fill = status)) +
  geom_boxplot(coef = 18) + 
  # guides(fill = FALSE) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(colour="black", size=12, face="italic")) +
  # scale_y_continuous(breaks = seq(0, 20000, by = 2500)) +
  coord_cartesian(ylim = c(0, 11.5)) +
  ylab("Abundance")


# # B: Dotplot per-sample variant
# ggplot(prop4, aes(x = gene, y = reads/1000, color = status)) +
#   geom_jitter() +
#   guides(fill = FALSE) +
#   theme_minimal() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x = element_text(colour="black", size=12, face="italic")) +
#   # scale_y_continuous(breaks = seq(0, 20000, by = 2500)) +
#   ylab("Abundance")


#################################################
#### GENE EXPRESSION IN EACH GROUP PER GENUS ####
#################################################

gene_reader_genus <- function(x,n1,n2,condition){
  # This function reads DIAMOND merged output from FG_reader.R
  x <- read.table(x)               
  x <- x[,-3]                                                                   # Remove strain
  colnames(x) <- c("sample","genus","gene_ID","reads","patient","status")
  
  x$status[x$status=="nonIBD"] <- "H"                                           # Replace "nonIBD" label
  x <- subset(x, x[ , 4] >= 5)                                                  # Remove genus/samples with less than 5 reads
  
  x <- x %>%
    select("genus","reads","status") %>%
    group_by(genus, status) %>% 
    summarise(abundance_genus = median(reads)) %>%                              # Average abundance in common condition / genus
    
    filter(status != "UC") %>%                                                  # Remove UC condition
    
    group_by(genus) %>%
    summarise(difference = diff(abundance_genus)) %>%                           # Substraction is (H - CD)
    
    mutate(direction = as.numeric(difference > 0)) %>%                          # Dummy variable to indicate the direction
    mutate(label_hjust = ifelse(difference < 0, 1.1, -0.1)) %>%                 # Dummy variable to indicate text position

    filter(difference > n1 | difference < n2) %>%                               # Keep most polar genus
                                                                                # Lower n1: more genus differentiated in H    
                                                                                # Lower n2: more genus differentiated in D
    {. ->> x}
  
  # Assign condition label to average differences
  x$direction[x$direction==1] <- "H"
  x$direction[x$direction==0] <- condition                                      # Disease status: e.g. "UC" or "CD"
  
  return(x)                        
  
}


diff_genus <- function(x,ylima=-120,ylimb=60,step,title){                            # Colors do not enter outside function (??)
  # https://stackoverflow.com/questions/42643189
  
  colors <- c("H" = "#00BA38", "CD" = "#F8766D")                                # https://stackoverflow.com/questions/8197559
  # #F8766D salmon, #619CFF cyan                                                # modify to cyan / salmon in these plots
  
  p <- ggplot(x, aes(x = reorder(genus, difference), y = difference,            # Sort values
                          fill = direction)) + 
    geom_bar(stat="identity", col = "black") +
    geom_text(aes(y = difference, label = genus, hjust = label_hjust)) +        # Add positioning information
    scale_fill_manual(values = colors) +                                        
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey80", 
                                            linetype = "dashed"),
          panel.grid.minor.x = element_blank()) +
    scale_y_continuous(breaks = seq(ylima, ylimb, step),
                       limits = c(ylima, ylimb)) +
    ylab("Abundance") +
    ggtitle(title) +
    coord_flip()
  
  return(p)
  
}

## TDCD (PROPIONATE KINASE) ANALYSIS ## 
tdcd_g <- gene_reader_genus(x="tdcd_genus.tsv",n1=15,n2=-10,condition="CD")
diff_genus(tdcd_g, ylima=-220, ylimb=180, step=40, title=expression(italic("tdcD")))
ggsave("../../results/MTG/tdcd.png",units="in", width=9, height=7.5, dpi=300)

## PDUW (PROPIONATE KINASE) ANALYSIS ## 
pduw_g <- gene_reader_genus(x="pduw_genus.tsv",n1=10,n2=-5,condition="CD")
diff_genus(pduw_g, ylima=-120, ylimb=260, step=40, title=expression(italic("pduW")))
ggsave("../../results/MTG/pduw.png",units="in", width=9, height=7.5, dpi=300)

## PCT (PROPIONATE COA TRF) ANALYSIS ## 
pct_g <- gene_reader_genus(x="pct_genus.tsv",n1=0.5,n2=-1,condition="UC")
diff_genus(pct_g, ylima=-30, ylimb=30,title=expression(italic("pct")))
ggsave("../../results/MTG/pct.png",units="in", width=8, height=7.5, dpi=300)

## PRPE (PROPIONATE COA LIGASE) ANALYSIS ## 
prpe_g <- gene_reader_genus(x="prpe_genus.tsv",n1=0.5,n2=-1,condition="UC")
diff_genus(prpe_g, ylima=-30, ylimb=30,title=expression(italic("prpE")))
ggsave("../../results/MTG/prpe.png",units="in", width=8, height=7.5, dpi=300)


# # Pairwise Wilcoxon test for every genus in the three conditions
# 
# # https://stackoverflow.com/questions/23042844/
# res <- lapply(split(tdcd_genus, tdcd_genus$genus), function(d){
#   pairwise.wilcox.test(d$reads, d$status, p.adjust.method = "BH")
# })
# 
# # Convert to DF and remove NA comparison (the one that compares itself)
# x <- as.data.frame(sapply(res, function(x) x$p.value))
# x <- x[complete.cases(x), ]
# 
# 
# # No rownames, test which comparison corresponds to which row and modify
# pairwise.wilcox.test(tdcd_genus[tdcd_genus$genus == "Anaerocolumna", 4], 
#                      tdcd_genus[tdcd_genus$genus == "Anaerocolumna", 6], 
#                      p.adjust.method = "BH")
# 
# rownames(x) <- c("HvsCD","CDvsUC","HvsUC")
# x <- x[-2,]                                                                     # remove UC-CD comparison 
# 
# 
# test <- lapply(split(tdcd_genus, tdcd_genus$genus), function(d){
#   median(d$reads)
# })
















##############################################
## EXPLORING SOURCES OF TECHNICAL VARIATION ##
##############################################

## 1. LIBRARY SIZE
samples <- read.table("metadata_MGX_reads.txt")[c(1,2,4,5)]
colnames(samples) <- c("sample","patient","status","kkreads")

ggplot(samples, aes(x=reorder(sample, kkreads),y=kkreads, fill=status)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  scale_y_continuous(n.breaks = 20) +
  geom_vline(xintercept = "HSMA33OT") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("samples") +
  ylab("million reads") +
  ggtitle("Library size of metagenomic samples")





