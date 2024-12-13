#Import cui et al dataset
#library(dplyr)
virus_data <- read.csv(here("data", "Cui_etal2014.csv"))
str(virus_data)# check number of columns and rows
#log transform the data
virus_data_log <- virus_data %>%
  mutate(Log_Genome_length = log(Genome.length..kb.),
         Log_Virion_volume = log(Virion.volume..nm.nm.nm.))

ggplot(virus_data_log, aes(x = Log_Genome_length, y = Log_Virion_volume)) +
  geom_point(color = "black", size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Regression line with confidence interval
  labs(
    title = "Log-Log Relationship Between Genome Length and Virion Volume",
    x = "Log(Genome Length kb)",
    y = "Log(Virion Volume(nm3)"
  ) +
  theme_minimal()

#Fit the linear model
Virus_plot <- lm(Log_Virion_volume ~  Log_Genome_length, data = virus_data_log)

# Model summary
summary(Virus_plot)

# Extract coefficients
intercept <- coef(Virus_plot)[1]  # ln(α)
slope <- coef(Virus_plot)[2]      # β

# transform intercept back from log to get alpha
exp(intercept)
slope


# Parameters from the model 
alpha <- 1181.807  
beta <- 1.515228   

# Genome length of 300 kb
L <- 300000  # Length in nucleotides

# Calculate the volume
V <- alpha * (L^beta)

# Print the result
cat("Estimated Volume of a 300 kb dsDNA Virus:", V, "nm³\n")

