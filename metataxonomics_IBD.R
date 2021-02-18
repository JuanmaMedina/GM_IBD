library("tidyverse")

setwd("../Desktop/IBBTEC/articulos/paper_IBD/methods/MTT/")

# Read metadata
samples <- read.table("metadata_MTT.txt")[c(1,2,4)]
colnames(samples) <- c("sample","patient","status")

# Add prefix to samples, as they are numeric and will need it later
samples$sample <- paste0('X', samples$sample)

# Replace "nonIBD" label
samples[samples=="nonIBD"] <- "H"

# N of each group
table(samples$status)


# TAXONOMIC ABUNDANCES
m <- read.csv("taxonomic_profiles.tsv", header = TRUE, sep = "\t")[,-1]         # First column are useless OTUs codes

# Parse taxonomy column
m$taxonomy <- str_remove_all(m$taxonomy, " __")

m %>%
  separate(col=taxonomy, 
           into=c("Domain","Phylum","Class","Order","Family","Genus"),          # Split on taxonomical columns
           sep = ";", remove=T) %>%  
  select(-c("Domain","Phylum","Class","Order","Family")) %>%                    # Keep only genus
  filter(!(Genus == "g" | Genus == "uncultured")) %>%                           # Remove taxonomies with undefined genus
  {. ->> m }

# Parsing crappy taxonomic annotations
m$Genus <- gsub('^_','',m$Genus)                                                # Remove weird initial underscores
m$Genus <- gsub("\\_.*","",m$Genus)                                             # Remove everything from underscore on
m <- m[- grep("ceae", m$Genus),]                                                # Remove families
m <- m[- grep("Family", m$Genus),]                                              # Remove additional undefined genus

# ## Parallel analysis to understand Price19 Bacteroides taxonomic information
# mb <- m[m$Genus == "Bacteroides", ]
# mb <- as_tibble(cbind(sample = names(mb), t(mb)))[-179,]
# mb <- merge(mb,samples[,c(1,3)],by = "sample")
# mb <- mb[,c(15,2:14)]
# mb[-1] <- lapply(mb[-1], as.numeric)
# mb <- mb %>% group_by(status) %>% summarise_all("median")
# mb <- mb[,c(1,2,3,7,8)]
# colnames(mb) <- c("status","speciesA","speciesB","speciesC","speciesD")
# mb <- gather(mb, species, abundance, speciesA:speciesD, factor_key=TRUE)
# ggplot(mb, aes(x = species, y = log(abundance+1,2), fill = status)) +
#   geom_bar(stat="identity", position=position_dodge()) + theme_minimal() + 
#   theme(axis.title.y=element_blank()) +
#   ylab("log2(Abundance)") + ggtitle("Bacteroides spp.") + coord_flip()

# Sum up abundance values of the same genus (first column as rownames)
m <- rowsum(x=m[,-179],group = m$Genus,)

# Transpose tibble object
m <- as_tibble(cbind(sample = names(m), t(m)))

# Merge both DFs by sample to get the disease status (NB: remove patients here)
m <- unique(merge(m,samples[,c(1,3)],by = "sample"))
m <- m[,c(329,2:328)]                                                           # status, N counts, remove samples

table(m$status)

m[-1] <- lapply(m[-1], as.numeric)                                          

# Average abundances in each sample and condition
m %>% 
  group_by(status) %>%
  summarise_all("median") %>%
  {. ->> m}

# # Log transform (this logs everything before substacting averages, more innacurate)
# m[,-1] <- log(m[,-1]+1, 10)

# Transpose, first row as colnames
m <- as_tibble(cbind(sample = names(m), t(m)))
names(m) <- as.character(unlist(m[1,])); m <- m[-1,]                            # First row as colnames
m <- m %>% remove_rownames %>% column_to_rownames(var="status")                 # First column as rownames


################################################################################


# Alistipes: 13, 123, 8 (CD, H, UC)

m %>%
  rownames_to_column("genus") %>%                                               # Preserve row names
  mutate_at(c("CD", "H", "UC"), as.numeric) %>%                                 # From character to numeric
  mutate(HminusCD = H - CD) %>%
  mutate(HminusUC = H - UC) %>%
  drop_na() %>%                                                                 # Weird NA generated
  filter(!(HminusCD==0 & HminusUC == 0)) %>%                                    # Remove genus with no quantification
  
  mutate(HminusCD = ifelse(HminusCD < 0, -log10(abs(HminusCD-1)), log10(HminusCD+1))) %>%    # log transform
  mutate(HminusUC = ifelse(HminusUC < 0, -log10(abs(HminusUC-1)), log10(HminusUC+1))) %>%
  
  mutate(directionCD = as.numeric(HminusCD > 0)) %>%                            # Dummy variables to indicate the direction
  mutate(directionUC = as.numeric(HminusUC > 0)) %>%
  
  mutate(label_hjust_CD = ifelse(HminusCD < 0, 1.1, -0.1)) %>%                  # Dummy variable to indicate text position
  mutate(label_hjust_UC = ifelse(HminusUC < 0, 1.1, -0.1)) %>% 
  
  mutate(directionCD = factor(ifelse(directionCD == 1, "H", "CD"))) %>%         # Conditional replacement 
  mutate(directionUC = factor(ifelse(directionUC == 1, "H", "UC"))) %>%         # https://stackoverflow.com/questions/35610437
  
  filter(HminusCD > 1.6 | HminusCD < -1) %>%                                       # Keep most polar genus (for logscale)
  # filter(HminusCD > 0.5 | HminusCD < 0) %>%                                   # Keep most polar genus (for normal scale)
  {. ->> mf}


ggplot(mf, aes(x = reorder(genus, HminusCD), y = HminusCD,                      # Sort values
              fill = directionCD)) +
  geom_bar(stat="identity", col = "black") +
  geom_text(aes(y = HminusCD, label = genus, 
                hjust = label_hjust_CD)) +                                      # Add positioning information
  scale_fill_manual(values = c("H" = "#00BA38", "CD" = "#F8766D")) +            # #F8766D salmon, #619CFF cyan
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80",
                                          linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(breaks = seq(-3.5, 5.5, 2),
                     limits = c(-3.5, 5.5)) +
  ggtitle(expression(italic("16S"))) +
  ylab("log10(Abundance)") +
  coord_flip()

ggplot(mf, aes(x = reorder(genus, HminusUC), y = HminusUC,                      # Sort values
               fill = directionUC)) +
  geom_bar(stat="identity", col = "black") +
  geom_text(aes(y = HminusUC, label = genus, 
                hjust = label_hjust_UC, size = 14)) +                           # Add positioning information
  scale_fill_manual(values = c("H" = "#00BA38", "UC" = "#619CFF")) +            # #F8766D salmon, #619CFF cyan
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80",
                                          linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(breaks = seq(-2.5, 4.5, 1),
                     limits = c(-2.5, 4.5)) +
  ylab("log10(Abundance)") +
  coord_flip()
