#########################################################################
###################### 1. LIBRARY OF PACKAGES ###########################
library(tidyverse)
library(ggplot2)
library(psych)
library(readr)
library(car)
library(ARTool)
library(vcd)
library(dplyr)
library(readxl)
library(janitor)
##########################################################################
###################### 2. CHARGE  DATA ################################

dades <- read_excel("C:/Users/pauss/Desktop/Apunts/TFG/Rstudio/excelR.xlsx")
view(dades)
str(dades)

#########################################################################
#######################3. Plots #########################################

dades <- dades %>%
  mutate(Lineage = factor(Lineage))

##### Gens x Llinatge ###
ggplot(dades, aes(x = Lineage, y = Gens, fill = Lineage)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Nombre de gens per llinatge",
    x = "Llinatge",
    y = "Nombre de gens"
  ) +
  theme(legend.position = "none")

#### Genoma x Llinatge #####


ggplot(dades, aes(x = Lineage, y = Genoma, fill = Lineage)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Mida del genoma per llinatge",
    x = "Llinatge",
    y = "Mida del genoma (bp)"
  ) +
  theme(legend.position = "none") 
+ scale_y_continuous(labels = scales::comma)

###### Prot x llinatge ####

ggplot(dades, aes(x = Lineage, y = Protein, fill = Lineage)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Nombre de proteïnes per llinatge",
    x = "Llinatge",
    y = "Nombre de proteïnes"
  ) +
  theme(legend.position = "none")

###### Resistencia pred. x llinatge #####


dades <- dades %>%
  rename(
    Num_antibiotics_resistents = `Resistència predita`
  ) %>%
  mutate(
    Num_antibiotics_resistents = as.numeric(Num_antibiotics_resistents),
    Lineage = factor(Lineage, levels = sort(unique(Lineage)))
  )

# Comprovació ràpida

dades <- dades %>%
  clean_names()

colnames(dades)


dades <- dades %>% 
  rename(num_antibiotics_resistents = resistencia_predita)%>%
  mutate(num_antibiotics_resistents = as.numeric(num_antibiotics_resistents),
    lineage = factor(lineage, levels = sort(unique(lineage)))
    
    
plot1 <- ggplot(dades, aes(
      x = lineage,
      y = num_antibiotics_resistents,
      fill = lineage
    )) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      theme_minimal() +
      labs(
        title = "Dispersió de la resistència predita per llinatge",
        x = "Llinatge",
        y = "Nombre d’antibiòtics resistents"
      ) +
      theme(legend.position = "none")
plot1
