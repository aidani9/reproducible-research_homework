# Reproducible research: version control and R

\## Question 4 a): Both walks begin at the origin (0, 0) and then diverge as they progress, each with its own unique trajectory despite running from the same code. The walk moves a fixed length in a random direction with each step, and the colour of the path changes with time. Since the walks were generated independently,  the 2 plots I made are very different from each other - with one circling around and ending up close to the origin, and the other straying very far from the origin. This illustrates the stochastic nature of these random walks and how if someone tried to use the same code they could generate very different results

4 b)  what is a random seed and how does it work: when we generate "random numbers" using a computer, they are not truly random- they give you a set sequence of "random" numbers depending an initial input (the seed). The numbers the computer gives you appear random but are actually deterministic based on the seed. The seed is like the starting point for the algorithm, if you use the same seed, you will get the same sequence of numbers. This is very important as it means others can use your code to get the same results that you did at the time.

4 c) Edited file in repo

4 d) ![edits to random walk script](https://github.com/user-attachments/assets/1a4ebe22-2ca3-4aa5-b9a4-a14912407bc1)

5 a) The data has 33 rows and 13 columns
5 b) Genome length and Virion volume span multiple orders of magnitude, so we can use a log transformaton to help the linear model capture these relationships better than the raw data, where the large calues may skew our results. 
#log transform the data
virus_data_log <- virus_data %>%
  mutate(Log_Genome_length = log(Genome.length..kb.),
         Log_Virion_volume = log(Virion.volume..nm.nm.nm.)) 

5 c) I found the exponent to be 1.515228 and scaling factor to be 1181.807 using the code below, and these results are the same as those in table 2 of the paper.  with the summary function I found the P value  2.279645e-10 for the exponent and  6.438498e-10 for the scaling factor, which are both far below our significance level of 0.05, so its very unlikely these results are due to chance. 
#Fit the linear model
> Virus_plot <- lm(Log_Virion_volume ~  Log_Genome_length, data = virus_data_log)
> # Extract coefficients
> intercept <- coef(Virus_plot)[1]  # ln(α)
> slope <- coef(Virus_plot)[2]      # β
> # transform intercept back from log to get alpha
> exp(intercept)
(Intercept) 
   1181.807 
> slope
Log_Genome_length 
         1.515228 

5 d) 
ggplot(virus_data_log, aes(x = Log_Genome_length, y = Log_Virion_volume)) +
  geom_point(color = "black", size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Regression line with confidence interval
  labs(
    title = "Log-Log Relationship Between Genome Length and Virion Volume",
    x = "Log(Genome Length kb)",
    y = "Log(Virion Volume(nm3)"
  ) +
  theme_minimal()

  5 e) using the code below, i found the estimated Volume of a 300 kb dsDNA Virus to be 235306564197 nm³ 
alpha <- 1181.807  
beta <- 1.515228   

# Genome length of 300 kb
L <- 300000  # Length in nucleotides

# Calculate the volume
V <- alpha * (L^beta)





## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points. First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers (also make sure that your username has been anonymised). All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   a) A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points) \
   b) Investigate the term **random seeds**. What is a random seed and how does it work? (5 points) \
   c) Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points) \
   d) Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points) 

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \alpha L^{\beta}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   a) Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)\
   b) What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points) \
   c) Find the exponent ($\beta$) and scaling factor ($\alpha$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points) \
   d) Write the code to reproduce the figure shown below. (10 points) 

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  e) What is the estimated volume of a 300 kb dsDNA virus? (4 points) 
