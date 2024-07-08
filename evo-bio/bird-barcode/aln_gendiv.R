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

## GOTTA FIGURE THIS OUT
library(forcats)
coi_data %>% arrange(pi) %>%
  mutate(binomial = factor(binomial, levels=binomial)) %>% 
  ggplot(aes(x=binomial,y=pi))+geom_point()+
  facet_wrap(~island)
