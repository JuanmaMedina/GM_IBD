library("tidyverse")

setwd("../Desktop/IBBTEC/articulos/paper_IBD/methods/MBX/")

# Read metadata
samples <- read.table("metadata_MBX.txt")[c(1,2,4)]
colnames(samples) <- c("sample","patient","status")

# Replace "nonIBD" label
samples[samples=="nonIBD"] <- "H"

# N of each group
table(samples$status)

# Extract propionate, lactate and butyrate
# tail -n+2 HMP2_metabolomics.csv | awk -F "," '{print $6}' | 
#     grep -i -n "formate\|acetate\|propionate\|butyrate\|lactate" 

# LINES AND METABOLITES THEY CONTAIN (### means not considered metabolite)
### 10:indoleacetate
# 84:phenyllactate
### 85:p-hydroxyphenylacetate
# 41363:2-aminobutyrate
# 41364:2-hydroxy-3-methylbutyrate
# 41390:butyrate
# 41424:imidazolelactate
# 41425:indole-3-propionate
# 41428:lactate
### 41444:phenylacetate
# 41448:propionate
# 57859:imidazole propionate

# Read target lines
mtb <- read.csv("HMP2_metabolomics.csv", header = TRUE, sep = ",")[c(84,41424,41428,     # lactate
                                                                     41363,41364,41390,  # butyrate
                                                                     41425,41448,57859   # propionate
                                                                     ), 
                                                                   -c(1,2,3,4,5,7)]

# Simpler variable and skip reading file in case of correction
m <- mtb
colnames(m)[1] <- "metabolite"

# There are some NAs in the DF. Consider as null values
m[is.na(m)] <- 0

# Test: these are ABSOLUTE values, as rows (metabolites) do not sum the same (NAs by 0 do not affect this)
rowSums(m[c(1,2,3),-1])
# Columns (samples) do not sum the same either (checked with several samples original file in bash)
# tail -n+2 HMP2_metabolomics.csv | awk -F "," '{s+=$100} END {print s}'

# Consider only the three target metabolites
m <- m[c(3,6,8),]

# First column as rownames
m %>% 
  remove_rownames %>% 
  column_to_rownames(var="metabolite") %>%
  {. ->> m}

# Transpose tibble object
m <- as_tibble(cbind(sample = names(m), t(m)))

# Merge both DFs by sample to get the disease status (NB: remove patients here)
m <- unique(merge(m,samples[,c(1,3)],by = "sample"))

# Convert values to numeric
m[, c(2,3,4)] <- sapply(m[, c(2,3,4)], as.numeric)

############################# 
########### PLOTS ###########
############################# 

### BOXPLOT FUNCTION ###
m_bp <- function(metabolite,a=0,b=2.5e06){
  
  x <- ggplot(m, aes(x=status, y=metabolite, fill=status)) + 
    geom_boxplot(factor=10) +
    # guides(fill = FALSE) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(colour="black", size=10)) +
    scale_y_continuous(limits = c(a,b), 
                       labels = function(x) format(x, scientific = TRUE)) +
    ylab("Abundance")
  
  return(x)
}

### BOXPLOTS ###
m_bp(m$butyrate,"Butyrate")
m_bp(m$propionate,a=0,b=2.5e06)
m_bp(m$lactate,"Lactate")

### VIOLINPLOT FUNCTION ###
m_vio <- function(metabolite,title){
  x <- ggplot(m, aes(x = status, y = metabolite, fill=status)) +
    geom_violin() + 
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(colour="black", size=10)) +
    ggtitle(title) +
    ylab("Abundance")
  
  return(x)
}

### VIOLINPLOTS ###
m_vio(m$butyrate,"Butyrate")
m_vio(m$propionate,"Propionate")
m_vio(m$lactate,"Lactate")

### SCATTERPLOT FUNCTION ###
m_sca <- function(metabolite,title){
  x <- ggplot(m, aes(x = status, y = metabolite, color=status)) +
    #geom_jitter(width = 0.4, height = 0, show.legend = FALSE) + 
    geom_jitter(show.legend = FALSE) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(colour="black", size=10)) +
    ggtitle(title) +
    ylab("Abundance")
  
  return(x)
}

### SCATTERPLOTS ###
m_sca(m$butyrate,"Butyrate")
m_sca(m$propionate,"Propionate")
m_sca(m$lactate,"Lactate")


m_den <- function(metabolite,a,b){                                              # Limits
  x <- ggplot(m, aes(x = metabolite)) +                                         # Re-scaling ( * 10^-6)                                       
    geom_density(aes(fill = status), alpha=0.7) + 
    theme_minimal() +
    guides(fill = FALSE) +
    ylab("Density") +
    scale_x_continuous(limits = c(a,b),
                         labels = function(x) format(x, scientific = TRUE)) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +

    xlab("Abundance")
  
  return(x)
}

### DENSITIES ###
m_den(m$butyrate,"Butyrate",b=20)
m_den(m$propionate,a=0,b=2.8e06)
m_den(m$lactate,"Lactate",b=3)


############################# 
######## STATISTICS #########
#############################

### KRUSKAL TEST AND DUNN TEST ###
# http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

# Is there any significant difference between the average metabolic measures in the
# three conditions? Use Kruskal-Wallis, non-parametric alternative to one-way ANOVA
# test which extends the two-samples Wilcoxon test to more than two groups. It is
# recommended when the assumptions of one-way ANOVA test are not met

kruskal.test(butyrate ~ status, data = m)   ## p-value = 0.1778
kruskal.test(propionate ~ status, data = m) ## p-value = 0.04036 <--
kruskal.test(lactate ~ status, data = m)    ## p-value = 0.8317

# As the p-value is less than the significance level 0.05, we can conclude that there are 
# significant differences in propionate levels between the three conditions. But we do not
# know which conditions are different

## Method 1: pairwise wilcoxon test (BH correction)
pairwise.wilcox.test(m$propionate, m$status, p.adjust.method = "BH")
aggregate(m[, 2:4], list(m$status), median)
aggregate(m[, 2:4], list(m$status), mean)

# The pairwise comparison shows that UC and healthy are significantly different 
# (p = 0.029 < 0.05)
ms <- m[m$propionate < 1.5e06,]
pairwise.wilcox.test(ms$propionate, ms$status, p.adjust.method = "BH")
aggregate(ms[, 2:4], list(ms$status), median)
aggregate(ms[, 2:4], list(ms$status), mean)

## Method 2: Dunn test
# https://cran.r-project.org/web/packages/dunn.test/dunn.test.pdf

library("dunn.test")

# If the Kruskal-Wallis test is significant, a posterior analysis can be performed to determine
# which levels of the independent variable (i.e. which conditions) differ from each other

# Dunn test is appropriate for groups with unequal numbers of observations, and can be used to 
# pin-point which specific means are significant from the others
# altp = T means we are working with p-vals (0.05 significance limit), not z-scores (0.025 
# significance limit). But they are equivalent

dunn.test(m$butyrate, m$status, method="bh", altp = T)
dunn.test(m$propionate, m$status, method="bh", altp = T)
dunn.test(m$lactate, m$status, method="bh", altp = T)

# The Dunn test shows that UC and healthy are significantly different in propionate production
# (p = 0.035 <= alpha == 0.05). 

mean(m[m$status == "H", "propionate"])
mean(m[m$status == "UC", "propionate"])
mean(m[m$status == "CD", "propionate"])

# Propionate production is significantly higher in H compared to UC


