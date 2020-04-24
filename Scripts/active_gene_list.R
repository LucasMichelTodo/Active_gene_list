library(readxl)
library(tidyverse)
library(reshape2)

## Read translation table
map <- read.csv('/Users/lucas/ISGlobal/PhD_Project/Data/oldnames_table.csv')
head(map)

## Import variantome table
variantome <- read_excel("~/ISGlobal/Variantome_RoviraGraells/3D7_Variantome_AllData_withGam.xls", sheet = 4)
head(variantome)
colnames(variantome)[1] <- "Old_id"
head(variantome)

df <- variantome %>%
  select(Old_id, AvgDif_12B = `1.2B`,
         AvgDif_10G = `10G`,
         AvgDif_3D7B = `3D7-B`,
         Red_12B = `Aver.2Higher1.2B.`,
         Red_10G = `Aver.2Higher10G.`,
         Red_3D7B = `Aver.2Higher3D7-B.`) %>%
  left_join(map, by='Old_id') %>%
  select(-Old_id) %>%
  group_by(Gene_id) %>% summarize_all(funs(mean))

head(df)
length(df$Gene_id)
length(unique(df$Gene_id))
df[!complete.cases(df),]

## Plot Red Signal Distribution
x <- df %>% select(Gene_id, starts_with("Red"))
mx <- melt(x)

p <- ggplot(mx, aes(x = value, color = variable, fill=variable)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 250, color = 'red')

ggsave('~/ISGlobal/PhD_Project/Plots/red_signal.png', p)

p <- ggplot(mx, aes(x = value, color = variable, fill=variable)) +
  geom_density(alpha = 0.2) +
  xlim(c(0,1000)) +
  geom_vline(xintercept = 250, color = 'red')

ggsave('~/ISGlobal/PhD_Project/Plots/red_signal_zoom.png', p)

## Clag Genes, for reference
df %>% filter(Gene_id == 'PF3D7_0302500')
df %>% filter(Gene_id == 'PF3D7_0302200')

## Plot Difference from Average distribution
x <- df %>% select(Gene_id, starts_with("Avg"))
mx <- melt(x)

p <- ggplot(mx, aes(x = value, color = variable, fill=variable)) +
  geom_density(alpha = 0.2)

print(p)

#ggsave('~/ISGlobal/PhD_Project/Plots/dif_from_avg.png', p)

p <- ggplot(mx, aes(x = value, color = variable, fill=variable)) +
  geom_density(alpha = 0.2) +
  xlim(c(-1.5,1.5))

ggsave('~/ISGlobal/PhD_Project/Plots/dif_from_avg_zoom.png', p)

## Set filtring parameters

thld_avgDif <- 0.4
thld_redSig <- 250

## Calculate Max dif from average

max_dif <- function(vect){
  mx <- max(vect, na.rm = T)
  mn <- min(vect, na.rm = T)
  if (is.infinite(mx) | is.infinite(mn)) {
    md <- NA
  } else {
    md <- mx - mn
  }
  return(md)
}

md <- df %>% select(starts_with('Avg')) %>%
  apply(1, max_dif)

df <- df
df['MaxAvgDif'] <- md

## Plot Avg MaxDif
head(df)

p <- ggplot(df, aes(x=MaxAvgDif)) +
  geom_density() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = 0.4, color = 'red')
ggsave('~/ISGlobal/PhD_Project/Plots/max_avg_dist.png', p)

## Remove rows in which AvgMaxdif == NA (means no avg for any strain) (79 rows)

table(is.na(df$MaxAvgDif))
df <- df %>% filter(!is.na(MaxAvgDif))
head(df)

## Set genes as active or inactively
## Set as inactive only if:
## - Red val. < thld
## - Differentialy expressed (max.avg.dist > thld) and dif from avg is negative

x <- df

x['Exp_12B'] <- 1
x[x$MaxAvgDif >= thld_avgDif & (!is.na(x$AvgDif_12B) & x$AvgDif_12B < 0),]$Exp_12B <- 0
x[x$Red_12B < thld_redSig,]$Exp_12B <- 0

x['Exp_10G'] <- 1
x[x$MaxAvgDif >= thld_avgDif & (!is.na(x$AvgDif_10G) & x$AvgDif_10G < 0),]$Exp_10G <- 0
x[x$Red_10G < thld_redSig,]$Exp_10G <- 0

x['Exp_3D7B'] <- 1
x[x$MaxAvgDif >= thld_avgDif & (!is.na(x$AvgDif_3D7B) & x$AvgDif_3D7B < 0),]$Exp_3D7B <- 0
x[x$Red_3D7B < thld_redSig,]$Exp_3D7B <- 0

write.csv(x, '~/ISGlobal/PhD_Project/active_genes.csv')
print(x, width=200)

## Clag Genes, for reference
print(x %>% filter(Gene_id == 'PF3D7_0302500'), width = 200)
print(x %>% filter(Gene_id == 'PF3D7_0302200'), width = 200)

print(x, width=200)
