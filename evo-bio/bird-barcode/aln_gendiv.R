## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
proj_path <- '~/Documents/GitHub/coding-resources/evo-bio/bird-barcode/'
setwd(proj_path)

## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
fasta <- list.files('data/bird/fasta/',full.names = T) %>% read_lines() %>% tibble() %>% rename(data = '.')


## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
code <- fasta %>% filter(str_detect(data,'>')) %>% mutate(data = str_replace(data,pattern = '>',replacement = '')) %>%
  rename(code = 'data')


## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
sequence <- fasta %>% filter(!str_detect(data,'>')) %>% rename(sequence = 'data')


## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
coi_data <- tibble(code,sequence) %>%
  mutate(code = str_replace(code,pattern='-[:digit:]$',replacement = '')) %>%
  mutate(code = str_replace(code,pattern='_[:digit:]$',replacement = '')) %>%
  distinct(code,.keep_all = TRUE)
  

## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
info <- read_csv('data/bird/bird_specimens.csv') %>% select(code,island,binomial)


## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
coi_data <- coi_data %>% left_join(info,by='code') %>% group_by(island,binomial) %>%
  summarize(sequence = list(sequence)) %>%
  group_by(island,binomial) %>% mutate(n = length(sequence[[1]])) %>% filter(n > 1)

## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
library(Biostrings)
coi_data <- coi_data %>% group_by(binomial) %>% mutate(sequence = list(DNAStringSet(unlist(sequence))))

## ----message=FALSE,warning=FALSE, results='hide'-------------------------------------------------------------------------------------------------------------------------------
library(muscle)
coi_data <- coi_data %>% mutate(aln = list(muscle::muscle(sequence[[1]])))

## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
library(pegas)
coi_data <- coi_data %>% mutate(pi = nuc.div(as.DNAbin(aln[[1]])))

## Plotting abundance and genetics

comm <- readRDS(file = 'data/bird/abundance/comm_counts_samples.rds')
comm_gen <- vector('list')
for (x in names(comm)) {
  comm_gen[[x]] <- comm[[x]] %>% 
    rename(sciName = 'binomial') %>%
    left_join(coi_data[which(coi_data$island==str_to_lower(x)),],by='binomial') %>% 
    select(!c(sequence,aln,island)) %>%
    mutate(binomial = fct_reorder(binomial,dplyr::desc(total)))
}
names(comm_gen) <- names(comm)

coeff <- 10^4 # Coefficient for second axis

i=1
ggplot(comm_gen[[i]],aes(x=binomial,y=total))+
  geom_bar(stat='identity',position = position_dodge(),fill='cadetblue3')+
  stat_smooth(aes(y=total, x=as.numeric(binomial)), method = lm,, formula = y ~ poly(x, 10),se = FALSE,
              col = 'red',size=0.4,linetype='dashed')+
  geom_point(aes(y=pi*coeff),color='purple',size=3.5)+
  scale_y_continuous(
    name = "Abundance",
    sec.axis = sec_axis(~./coeff)
  )+
  labs(x='Species')+
  theme_minimal()+
  #theme_light()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=17),
        axis.title.x = element_blank(),
        #axis.title.x = element_text(size=18, vjust = 0.1),
        axis.title.y = element_text(size=20, vjust = 1.8),
        plot.title = element_text(size=26))+
  ggtitle(names(comm_gen)[i])
ggsave(paste0('output/',names(comm_gen)[i],'_SAD_pi.png'),width = 9,height = 7,dpi=1200)

## Hill numbers

x <- lapply(comm,select,c(1,2))
# Getting spp-site matrix
getMat <- function (x) {
  sppmat <- x[[1]]
  colnames(sppmat)[2] <- names(x)[1]
  if (length(x)>1) {
    for (i in 2:length(x)) {
      colnames(x[[i]])[2] <- names(x)[i]
      sppmat <- sppmat %>% full_join(x[[i]])
    }
  }
  sppmat <- sppmat %>% mutate(across(2:ncol(sppmat),~replace_na(.x,0)))
  allspp <- sppmat$sciName
  sppmat <- data.frame(sppmat[,2:ncol(sppmat)])
  rownames(sppmat) <- allspp
  return(t(sppmat))
}
sppmat <- getMat(x)
require(hillR)
taxahill <- sapply(0:3,hill_taxa, comm = sppmat)
colnames(taxahill) <- paste0('abund_h',0:3)
