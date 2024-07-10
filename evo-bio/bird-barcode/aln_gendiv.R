## ----message=FALSE,warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
proj_path <- '/Volumes/G-DRIVE USB/Rilquer/github_repos/coding-resources/evo-bio/bird-barcode/'
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
comm <- lapply(1:length(comm),function(x){return(mutate(comm[[x]],island = rep(names(comm)[x],nrow(comm[[x]]))))}) %>%
  do.call(what = rbind.data.frame) %>% mutate(island = str_to_lower(island)) %>%
  rename(sciName = 'binomial') %>%
  left_join(coi_data,by=c('island','binomial')) %>%
  select(!c(sequence,aln,)) %>%
  #replace_na(list(n = 0, pi = 0)) %>%
  mutate(binomial = fct_reorder(binomial,dplyr::desc(total)))
coeff <- 10^3.2 # Coefficient for second axis
library(ggthemes)

comm %>%
  ggplot(aes(x=binomial,y=total,group=island))+
  geom_bar(stat='identity',position = position_dodge(),fill='cadetblue3')+
  stat_smooth(aes(y=total, x=binomial), method = lm, formula = y ~ poly(x, 10), se = FALSE,
              col = 'lightblue',size=0.5)+
  geom_point(aes(y=pi*coeff),color='purple')+
  facet_wrap(~island)+
  scale_y_continuous(
    name = "Abundance",
    sec.axis = sec_axis(~./coeff, name="Nuc. Div.")
  )+
  labs(x='Species')+
  #theme_hc()
  theme_few()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())
