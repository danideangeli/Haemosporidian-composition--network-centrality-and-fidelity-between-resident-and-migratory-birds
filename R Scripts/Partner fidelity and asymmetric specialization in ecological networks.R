title: "Partner fidelity and asymmetric specialization in ecological networks"
subtitle: "Supplementary material"
#author: "Matthew A. Barbour"
#date: "`r Sys.Date()`"
output:
  pdf_document:
  highlight: haddock
keep_tex: true
toc: true
number_sections: false
df_print: paged
fig_caption: true
citation_package: natbib
includes:  
  in_header: Am-Nat-preamble-latex.tex
biblio-style: amnat
bibliography: references
fontsize: 11pt
fontfamily: mathpazo
documentclass: article
linkcolor: black
urlcolor: black
citecolor: black
---
  
  ```{r setup, include=FALSE}
set.seed(34) # make sure simulations give same output
# load required libraries
library(tufte)
library(tidyverse)
library(brms)
library(cowplot)
library(broom)
library(kableExtra)
library(gmodels)        # make custom contrasts
library(codingMatrices) # visualize contrast matrix
library(lme4)           # for glmer models in permutation routine
# invalidate cache when the tufte version changes
knitr::opts_chunk$set(tidy = FALSE, cache.extra = packageVersion('tufte'))
options(htmltools.dir.version = FALSE)
# load and manage data
full.df <- read_csv("data/dataset.csv") %>%
  rename(resource_sp = species_1,
         consumer_sp = species_2,
         r_ND = normalized_degree_sp_1,
         c_ND = normalized_degree_sp_2) %>%
  mutate(network_id = as.factor(network_id),
         id_pair = as.factor(id_pair),
         type = as.factor(type),
         subtype = as.factor(subtype)) 
weighted_subset.df <- read_csv("data/weighted_subset.csv") %>%
  rename(resource_sp = species_1,
         consumer_sp = species_2,
         r_ND = normalized_degree_sp_1,
         c_ND = normalized_degree_sp_2) %>%
  mutate(network_id = as.factor(network_id),
         id_pair = as.factor(id_pair),
         type = as.factor(type),
         subtype = as.factor(subtype), 
         log.sum_r = log(sum_1),
         log.sum_c = log(sum_2))
# set the number of cores for model
options(mc.cores=parallel::detectCores ()) 
```

\newpage

All data and code to reproduce the analyses presented below are available on GitHub (https://github.com/mabarbour/partner_fidelity) and have been archived on Zenodo (https://zenodo.org/badge/latestdoi/107263795).

## Dataset

Here's a summary of the full dataset we used in our analyses (table \ref{tab:full-dataset}).

```{r full-dataset, echo=FALSE}
sample_n(full.df, size=10) %>%
  select(connected, type, subtype, r_ND, c_ND, resource_sp, consumer_sp, id_pair, network_id) %>%  
  mutate(r_ND = round(r_ND, 2), 
         c_ND = round(c_ND, 2)) %>%
  knitr::kable(., caption = "Random sample of 10 rows from the full dataset.", booktabs = T) %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))
```

- **connected**: observed species interaction (0 = no interaction observed, 1 = interaction observed)
- **type**: type of interaction (A = antagonistic; M = mutualistic)
- **subtype**: specific type of interaction (`r unique(full.df$subtype)`)
- **r_ND**: normalized degree of resource species (mean = `r round(mean(full.df$r_ND), 2)`, SD = `r round(sd(full.df$r_ND), 2)`)
- **c_ND**: normalized degree of consumer species (mean = `r round(mean(full.df$c_ND), 2)`, SD = `r round(sd(full.df$c_ND), 2)`)
- **resource_sp**: resource species ID (n = `r length(unique(full.df$resource_sp))`)
- **consumer_sp**: consumer species ID (n = `r length(unique(full.df$consumer_sp))`)
- **id_pair**: consumer-resource interaction ID (n = `r length(unique(full.df$id_pair))`)
- **network_id**: ecological network ID (n = `r length(unique(full.df$network_id))`; labels match those on the [Web of Life](http://www.web-of-life.es/))


We also analyze a subset of these data that only consists of weighted networks (table \ref{tab:weighted-subset-dataset}), which we use to investigate the effects of sampling effort. 

```{r weighted-subset-dataset, echo=FALSE}
sample_n(weighted_subset.df, size=10) %>%
  select(connected, type, subtype, r_ND, c_ND, sum_shared, resource_sp, consumer_sp, id_pair, network_id) %>%  
  mutate(r_ND = round(r_ND, 2), 
         c_ND = round(c_ND, 2)) %>%
  kable(., caption = "Random sample of 10 rows from the subset of weighted networks.", booktabs = T) %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))
```

- **connected**: 0 = no interaction observed, 1 = interaction observed
- **type**: A = antagonistic, M = mutualistic
- **subtype**: `r unique(weighted_subset.df$subtype)`
- **r_ND**: mean = `r round(mean(weighted_subset.df$r_ND), 2)`, SD = `r round(sd(weighted_subset.df$r_ND), 1)`
- **c_ND**: mean = `r round(mean(weighted_subset.df$c_ND), 2)`, SD = `r round(sd(weighted_subset.df$c_ND), 1)`
- **sum_shared**: mean = `r round(mean(weighted_subset.df$log.sum_r), 1)`, SD = `r round(sd(weighted_subset.df$log.sum_r), 1)`
- **resource_sp**: n = `r length(unique(weighted_subset.df$resource_sp))`
- **consumer_sp**: n = `r length(unique(weighted_subset.df$consumer_sp))`
- **id_pair**: n = `r length(unique(weighted_subset.df$id_pair))`
- **network_id**: n = `r length(unique(weighted_subset.df$network_id))`


Note that the general structure of the subset of data is the same as the full dataset, except we have added one new variable. This variable corresponds to the number of observed interactions (sum of interaction frequency, **sum_shared**) for a consumer-resource pair in a network. This variable gives information on sampling effort<!--, which we later use to test whether they can explain our results-->. Note also that we only have three interaction **subtype**s now (no weighted networks for PlantHerbivore).

Prior to our analyses, we scaled normalized degree (mean = 0, SD = 1) for both consumers and resources. This allowed us to interpret both main effects and statistical interactions within the same model, and also allowed us to easily compare their effect sizes [@Schielzeth2010]. We also scaled the logarithm of observed interactions for both consumers and resources.

```{r scale variables}
full.df <- full.df %>%
  mutate(sc.r_ND = scale(r_ND),
         sc.c_ND = scale(c_ND))
weighted_subset.df <- weighted_subset.df %>%
  mutate(sc.r_ND = scale(r_ND),
         sc.c_ND = scale(c_ND),
         sc.log.sum_r = scale(log.sum_r),
         sc.log.sum_c = scale(log.sum_c))
```


For **type**, we created a contrast so that the intercept term represents the average probability (Ave) of an interaction across mutualistic (M) and antagonistic (A) interactions, and this coefficient represents the effect of mutualistic (relative to antagonistic) interactions (table \ref{tab:type-contrasts}).

```{r type-contrasts, echo=F}
# create contrast matrix
type_matrix <- rbind("_M.vs.A" = c(-1, 1))
# make contrasts for R
type_contrast <- make.contrasts(type_matrix)
# clarify description of contrasts
show_type_contrast <- mean_contrasts(type_contrast)
# add column names for clarity
attr(show_type_contrast, "dimnames")[[2]] <- levels(full.df$type)
# show contrasts
show_type_contrast %>%
  kable(., caption = "Contrasts for \\textbf{type}.", booktabs = T) %>%
  kable_styling(latex_options = c("hold_position"))
# set contrasts for both datasets
contrasts(full.df$type) <- type_contrast
contrasts(weighted_subset.df$type) <- type_contrast
```


```{r type contrasts for full, eval=FALSE, include=FALSE}
# create contrast matrix
type_matrix <- rbind("M.vs.A" = c(-1/2, -1/2, 1/2, 1/2),
                     "Poll.vs.Disp" = c(0, 0, 1, -1),
                     "Herb.vs.Para" = c(-1, 1, 0, 0))
# make contrasts for R
type_contrast <- make.contrasts(type_matrix)
# clarify description of contrasts
show_type_contrast <- mean_contrasts(type_contrast)
# add column names for clarity
attr(show_type_contrast, "dimnames")[[2]] <- levels(full.df$type)
# show contrasts
show_type_contrast
# set contrasts
contrasts(full.df$type) <- type_contrast
```


To test for the effect of network **subtype** within each **type**, we create two new variables:

```{r}
full.df <- full.df %>%
  mutate(subtype_Herb.vs.Para = ifelse(type == "M", 0,
                            ifelse(subtype == "PlantHerbivore", 1/2, -1/2)),
         subtype_Poll.vs.Disp = ifelse(type == "A", 0,
                            ifelse(subtype == "PlantPollinator", 1/2, -1/2)))
```

```{r eval=F, include=F}
summary(lm(connected ~ type + subtype_Herb.vs.Para + subtype_Poll.vs.Disp, data=full.df))
group_by(full.df, type, subtype) %>%
  summarise_at(vars(connected), list(mean))
mean(c(0.529, 0.591, 0.620, 0.744))
mean(c(0.620, 0.744)) - mean(c(0.529, 0.591))
0.591 - 0.529
0.62 - 0.744
```

Now, **subtype_Herb.vs.Para** tests for an effect of herbivory (relative to parasitism), and **subtype_Poll.vs.Disp** tests for an effect of pollination (relative to seed dispersal).

Since we do not have multiple **subtype**s for antagonistic interactions in the data subset, we only fit a **subtype** contrast for mutualistic interactions.

```{r type contrasts for subset, eval=F, include=F}
# create contrast matrix
subset.type_matrix <- rbind("M.vs.A" = c(-1, 1/2, 1/2),
                            "Poll.vs.Disp" = c(0, 1, -1))
                          # "Herb.vs.Para" = c(-1, 1, 0, 0)
# make contrasts for R
subset.type_contrast <- make.contrasts(subset.type_matrix)
# clarify description of contrasts
show_subset.type_contrast <- mean_contrasts(subset.type_contrast)
# add column names for clarity
attr(show_subset.type_contrast, "dimnames")[[2]] <- levels(weighted_subset.df$type) # should be subtype?
# show contrasts
show_subset.type_contrast
# set contrasts
contrasts(weighted_subset.df$type) <- subset.type_contrast
```

```{r}
weighted_subset.df <- weighted_subset.df %>%
  mutate(subtype_Poll.vs.Disp = ifelse(type == "A", 0,
                            ifelse(subtype == "PlantPollinator", 1/2, -1/2)))
```

\ 

## Statistical Models

We fit the following statistical model to our full dataset:

```{r full model formula}
interaction.formula <- brmsformula(
  connected ~ type + subtype_Herb.vs.Para + subtype_Poll.vs.Disp + sc.r_ND + sc.c_ND + 
    type:sc.r_ND + type:sc.c_ND + sc.r_ND:sc.c_ND +
    type:sc.r_ND:sc.c_ND +
    (1 | resource_sp) + (1 | consumer_sp) + (1 | id_pair) + (1 | network_id),
  family = bernoulli(link = "logit")
)
```

\ 

and a similar model to the data subset:

```{r subset model formula}
interaction_subset.formula <- brmsformula(
  connected ~ type + subtype_Poll.vs.Disp + sc.r_ND + sc.c_ND + 
    type:sc.r_ND + type:sc.c_ND + sc.r_ND:sc.c_ND +
    type:sc.r_ND:sc.c_ND +
    (1 | resource_sp) + (1 | consumer_sp) + (1 | id_pair) + (1 | network_id),
  family = bernoulli(link = "logit") 
) # note that we do not include subtype_Herb.vs.Para, since we do not have
# PlantHerbivore interactions in the weighted subset.
```

These models analyze the probability of a species interaction as a function of all main effects, as well as two- and three-way interactions between the type of interaction (**type**) and scaled normalized degree of resources (**sc.r_ND**) and consumers (**sc.c_ND**)(fixed effects). Note that we only test for main effects of network subtypes within the type of interaction. This model also allows the probability of interactions to vary among resource and consumer species, the unique consumer-resource pair, as well as among unique ecological networks (random effects). Note that by including the unique consumer-resource pair as a random effect, we are testing how the different factors in our model influence partner fidelity (i.e., probability of two species to interact when they co-occur). Note also that specifying a Bernoulli distribution on the logit scale for the error distribution in our model makes it so that the fixed and random effects in our statistical model represent effects on the log odds of partner fidelity.


\ 

## Choosing Priors

Given the number of random effects in our model, we used a Bayesian approach to estimate our model's parameters (as advised by [Bolker et al. 2008](https://www.sciencedirect.com/science/article/pii/S0169534709000196)). This required us to choose prior distributions for the fixed and random effects in our model.

```{r inspect priors, include=FALSE}
# see which priors we need to set.
get_prior(interaction.formula, data = full.df, family = bernoulli()) %>%
  select(prior:group) %>%
  kable(., caption = "Default priors for full model.", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))
```

We specified the following priors in our models:
  
  ```{r set priors}
interaction.priors <- c(
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A"),
  set_prior("normal(0,2)", class = "b", coef = "subtype_Herb.vs.Para"),
  set_prior("normal(0,2)", class = "b", coef = "subtype_Poll.vs.Disp"),
  set_prior("normal(1,2)", class = "b", coef = "sc.r_ND"),
  set_prior("normal(1,2)", class = "b", coef = "sc.c_ND"),
  set_prior("normal(0,2)", class = "b", coef = "sc.r_ND:sc.c_ND"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.r_ND"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.c_ND"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.r_ND:sc.c_ND"),
  set_prior("normal(0,2)", class = "sd"))
interaction_subset.priors <- c(
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A"),
  set_prior("normal(0,2)", class = "b", coef = "subtype_Poll.vs.Disp"),
  set_prior("normal(1,2)", class = "b", coef = "sc.r_ND"),
  set_prior("normal(1,2)", class = "b", coef = "sc.c_ND"),
  set_prior("normal(0,2)", class = "b", coef = "sc.r_ND:sc.c_ND"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.r_ND"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.c_ND"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.r_ND:sc.c_ND"),
  set_prior("normal(0,2)", class = "sd"))
sampling_effect.priors <- c(
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A"),
  set_prior("normal(0,2)", class = "b", coef = "subtype_Poll.vs.Disp"),
  set_prior("normal(1,2)", class = "b", coef = "sc.log.sum_r"),
  set_prior("normal(1,2)", class = "b", coef = "sc.log.sum_c"),
  set_prior("normal(0,2)", class = "b", coef = "sc.log.sum_r:sc.log.sum_c"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.log.sum_r"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.log.sum_c"),
  set_prior("normal(0,2)", class = "b", coef = "type_M.vs.A:sc.log.sum_r:sc.log.sum_c"),
  set_prior("normal(0,2)", class = "sd"))
```

Choosing priors requires thinking on the scale of the model. Therefore, we first give a description of what each term quantifies in these models before justifying our choices. Note that all of these fixed effects (except for **Intercept**) can be interpreted as an effect on the log-odds of partner fidelity:
  
  - **Intercept** = log-odds of partner fidelity across interaction types at the mean normalized degree of resources and consumers.
- **type_M.vs.A** = effect of a mutualistic interaction (relative to antagonistic).
- **subtype_Herb.vs.Para** = effect of herbivory relative to parasitism.
- **subtype_Poll.vs.Disp** = effect of pollination relative to seed disperser. 
- **sc.r_ND** = effect of 1 SD increase in resource normalized degree across interaction types.
- **sc.c_ND** = effect of 1 SD increase in consumer normalized degree across interaction types.
- **type_M.vs.A:sc.r_ND** = effect of a mutualistic interaction (relative to antagonistic) on a 1 SD increase in resource normalized degree.
- **type_M.vs.A:sc.c_ND** = effect of a mutualistic interaction (relative to antagonistic) on a 1 SD increase in consumer normalized degree.
- **sc.r_ND:sc.c_ND** = nonadditive effect of 1 SD increase in resource and consumer normalized degree across interaction types. A positive value indicates that symmetry in partners normalized degree increases partner fidelity, whereas a negative value indicates that *asymmetry* in partners normalized degree increases partner fidelity (see justification for this in section [**Asymmetry effect**](#asymmetry-effect)).
  - **type_M.vs.A:sc.r_ND:sc.c_ND** = nonadditive effect of a mutualistic interaction (relative to antagonistic) on a 1 SD increase in resource and consumer normalized degree. A positive value indicates that mutualistic interactions enhance the effects of symmetry on partner fidelity, whereas a negative value indicates that mutualistic interactions enhance the effects of *asymmetry*.
  
  For each of our random effects (**resource_sp**, **consumer_sp**, **id_pair**, and **network_id**), our model estimates the SD that describes the assumed normal distribution in their effect sizes.
  
  Below, we give a justification for each of the prior distributions we chose for these parameters.
  
  ### Type of interaction: **type_M.vs.A**
  
  Let's take a hypothetical example, where the effect of **type_M.vs.A** on partner fidelity is really large, ranging from a probability of 0.25 to 0.75. The example below shows that this corresponds to a logistic regression coefficient of ~2.2.

```{r type_M.vs.Acoefficient}
# probability of a mutualistic interaction when two species co-occur
pM <- 0.75  
# probability of an antagonistic interaction when two species co-occur
pA <- 0.25 
# logistic regression coefficient, which is the difference in log-odds 
# for 1 unit increase in the predictor
(mu_mutualist <- log(pM / (1-pM)) - log(pA / (1-pA))) 
```


We would expect such a large effect to be unlikely though. Also, we have no strong prior expectation as to whether mutualistic or antagonistic interactions will have a positive or negative effect on partner fidelity. So let's explore what a normal distribution looks like when the coefficient is centered on zero, but the variance is large enough to allow for a large effect (**type_M.vs.A** = 2.2, denoted by dotted line below) if there is enough evidence to support it.
  
  ```{r simulate prior distribution for type, echo=FALSE, fig.cap="Possible prior distributions for the effect of interaction type."}
  samp <- 1000
  series <- seq(-10, 10, length.out = samp)
  sd_type <- data.frame(series = rep(series, 4),
                        density = c(dnorm(series, 0, 0.5),
                                    dnorm(series, 0, 1),
                                    dnorm(series, 0, 2),
                                    dnorm(series, 0, 4)),
                        SD = c(rep("0.5", samp),
                               rep("1", samp),
                               rep("2", samp),
                               rep("4", samp)))
  ggplot(sd_type, aes(x = series, y = density, color = SD)) + 
    geom_line() + 
    xlab("Coefficient estimate") +
    ylab("Density") +
    scale_fill_discrete(name = "SD", labels = c(0.5, 1, 2, 4)) +
    geom_vline(xintercept = mu_mutualist, linetype = "dotted") # big effect
  ```
  
  Based on these distributions, we think a normal prior with mean=0 and SD=2 would be appropriate. This creates a regularizing prior where the mass of the distribution is centered on zero (i.e., no effect), but allows for the model to detect large effects if the data supports it.
  
  ### Subtype of interaction: **subtype_Herb.vs.Para** and **subtype_Poll.vs.Disp**
  
  We expect similar effects for different subtypes of interaction as we would for different interaction types. Therefore, we choose to use the same prior: normal(mean=0, SD=2).
  
  ### Scaled normalized degree: **sc.r_ND** and **sc.c_ND**
  
  <!--We expect for the normalized degree of a species to always have a positive relationship with the probability of an interaction. Biologically, it does not make sense that this relationship could ever be negative, which is why we will impose a lower bound of zero on this parameter. By default, we would expect this relationship to be 1:1 between the probability of observing an interaction and the normalized degree of a species. -->Since the normalized degree of a species defines its probability of interacting with a co-occuring partner, regardless of its identity, our prior expectation is that there is a 1:1 relationship between the normalized degree of a species and its probability of interacting with another species. In other words, the probability of observing an interaction is the same as the normalized degree. 
  
  Let's get some intuition as to what this prior looks like when normalized degree is on a standardized scale (1 SD above the mean) and the response is in terms of log odds, which matches the assumptions of our model. Note that we chose to only look 1 SD above the mean because this reflects how the coefficient is estimated in the model (effect of 1 unit increase in predictor variable). We do this first for resources: 

```{r resource ND prior}
# probability of average species interacting
pMean_r <- mean(full.df$r_ND) 
# probability of a relatively generalized species interacting
pGen_r <- mean(full.df$r_ND) + sd(full.df$r_ND) 
# logistic regression coefficient, which is the difference in log-odds 
# for a 1 SD unit increase in a species normalized degree.
(mu_sc.r_ND <- log(pGen_r / (1-pGen_r)) - log(pMean_r / (1-pMean_r))) 
```

and then for consumers: 

```{r consumer ND prior}
## Consumer species 
# probability of average species interacting
pMean_c <- mean(full.df$c_ND) 
# probability of a relatively generalized species interacting
pGen_c <- mean(full.df$c_ND) + sd(full.df$c_ND) 
# logistic regression coefficient, which is the difference in log-odds 
# for a 1 SD increase in a species normalized degree.
(mu_sc.c_ND <- log(pGen_c / (1-pGen_c)) - log(pMean_c / (1-pMean_c))) 
```

\ 

These examples indicate that we should have a prior expectation for the mean estimate of these coefficients to be ~1. How much variance around this coefficient might we expect though? Let's simulate a normal distribution with mean=1 (dotted line below), but different variances:
    
    ```{r simulate prior distribution for normalized degree, echo=FALSE, fig.cap="Potential prior distributions for the effect of scaled normalized degree."}
  samp <- 1000
  sd_sp1 <- data.frame(series = series,
                       density = c(dnorm(series, mu_sc.r_ND, 0.5),
                                   dnorm(series, mu_sc.r_ND, 1),
                                   dnorm(series, mu_sc.r_ND, 2),
                                   dnorm(series, mu_sc.r_ND, 4)),
                       SD = c(rep("0.5", samp),
                              rep("1", samp),
                              rep("2", samp),
                              rep("4", samp)))
  ggplot(sd_sp1, aes(x = series, y = density, color = SD)) + 
    geom_line() + 
    xlab("Coefficient estimate") +
    ylab("Density") +
    scale_fill_discrete(name = "SD", labels = c(0.5, 1, 2, 4)) +
    geom_vline(xintercept = 1, linetype = "dotted")
  ```
  
  Most standardized logistic regression coefficients are less than 5 [@Gelman2008], so we feel that specifying a standard deviation of 2 is reasonable since it more than covers the range from 0 to 5. It will also make the data "work" for larger values, since the mass of the distribution is centered on 1.
  
  ### Statistical interactions between interaction type and normalized degree: **type_M.vs.A:sc.r_ND**, **type_M.vs.A:sc.c_ND**
  
  We specified the same regularizing prior as for **type_M.vs.A** (i.e., normal(mean=0, SD=2)), because we have no strong prior expectation for whether interaction type will have a positive or negative effect on these relationships, and we want to make the data "work" for any strong effects.
  
  ### Partner (a)symmetry effect: sc.r_ND:sc.c_ND {#asymmetry-effect}
  
  ```{r eval=F, include=F}
  ND.x.ND_interpretation <- function(b_sc.r_ND = 1, b_sc.c_ND = 1, b_sc.r_ND.x.sc.c_ND, delta.r_ND, delta.c_ND){
    
    log.odds.connected <- b_sc.r_ND * delta.r_ND + 
      b_sc.c_ND * delta.c_ND + 
      b_sc.r_ND.x.sc.c_ND * delta.r_ND * delta.c_ND
    
    return(round(log.odds.connected, 2))
  }
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = 0.5, delta.r_ND = 0.5, delta.c_ND = 0.5)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = 0.5, delta.r_ND = 0.75, delta.c_ND = 0.25)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = 0.5, delta.r_ND = -0.5, delta.c_ND = -0.5)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = 0.5, delta.r_ND = -0.75, delta.c_ND = -0.25)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = -0.5, delta.r_ND = 0.5, delta.c_ND = 0.5)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = -0.5, delta.r_ND = 0.75, delta.c_ND = 0.25)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = -0.5, delta.r_ND = -0.5, delta.c_ND = -0.5)
  ND.x.ND_interpretation(b_sc.r_ND.x.sc.c_ND = -0.5, delta.r_ND = -0.75, delta.c_ND = -0.25)
  summary(full.df$sc.r_ND)
  summary(full.df$sc.c_ND)
  # covers >85% of the data
  coda::HPDinterval(as.mcmc(full.df$sc.r_ND), prob =0.85)
  coda::HPDinterval(as.mcmc(full.df$sc.c_ND), prob =0.85)
  ```
  
  Although it is not intuitive, the sign of the statistical interaction between consumer and resource normalized degree indicates whether partner symmetry or asymmetry increases partner fidelity. To illustrate this, let's first assume that the main effects of scaled normalized degree for both resources and consumers equals 1 (also our prior expectation). For simplicity, let's assume that the statistical interaction is either -1 or +1. With these coefficients, we can now calculate how a change in scaled normalized degree for both resources and consumers affects the log-odds of an interaction. Below, we explore these effects $\pm$ 1 SD for the scaled normalized degree of both resources and consumers.  
  
  ```{r get data for NDxND interaction, echo=T}
  sim.log.odds <- expand.grid(
    b_sc.r_ND = 1, 
    b_sc.c_ND = 1,
    `b_sc.r_ND:sc.c_ND` = c(-1, 1),
    delta.sc.r_ND = seq(-1, 1, 0.01),
    delta.sc.c_ND = seq(-1, 1, 0.01)) %>%
    mutate(log.odds.connected = 
             b_sc.r_ND * delta.sc.r_ND +                          # main effect of sc.r_ND
             b_sc.c_ND * delta.sc.c_ND +                          # main effect of sc.c_ND
             `b_sc.r_ND:sc.c_ND` * delta.sc.r_ND * delta.sc.c_ND, # statistical interaction
           main.log.odds = b_sc.r_ND * delta.sc.r_ND + b_sc.c_ND * delta.sc.c_ND,
           diff.log.odds = log.odds.connected - main.log.odds)
  ```
  
  ```{r plot-NDxND, echo=F, fig.cap="\\label{fig:plot-NDxND}Plot of how the sign of the statistical interaction between resource and consumer normalized degrees (b_sc.r_ND:sc.c_ND) determines how symmetry (right panel) or asymmetry (left panel) in partners normalized degree affects partner fidelity.", warning=FALSE}
  labels_2D <- list("-1"="Coefficient of normalized degree interaction\n(b_sc.r_ND:sc.c_ND) = -1", "1"="Coefficient of normalized degree interaction\n(b_sc.r_ND:sc.c_ND) = 1")
  labeller_2D <- function(variable, value) {
    return(labels_2D[value])
  }
  ggplot(sim.log.odds, aes(x = delta.sc.r_ND, y = delta.sc.c_ND, z=diff.log.odds)) +
    geom_raster(aes(fill = diff.log.odds)) +
    facet_wrap(~factor(`b_sc.r_ND:sc.c_ND`), labeller = labeller_2D) +
    scale_fill_viridis_c(name = expression(paste(Delta, "log odds"))) +
    xlab(expression(paste("Effect of resource normalized degree (", Delta, "sc.r_ND)"))) +
    ylab(expression(paste("Effect of consumer normalized degree (", Delta, "sc.c_ND)"))) +
    theme_cowplot(font_size = 10)
  ```
  
  By inspecting figure \ref{fig:plot-NDxND}, one can see how the sign of the statistical interaction between resource and consumer normalized degree affects the log odds of an interaction. When the sign is positive (right panel), symmetric normalized degrees increase the probability of an interaction (relative to main effects). In contrast, when the sign is negative (left panel), asymmetric normalized degrees increase the probability of an interaction. For example, if a relatively specialized consumer ($\Delta$sc.c_ND = -1) interacts with a relatively generalist resource ($\Delta$sc.r_ND = +1), there is an increase in the log-odds of an interaction (yellow spot in lower right corner of left panel in fig. \ref{fig:plot-NDxND}).
  
  Since we had no prior expectation as to how symmetry or asymmetry would affect partner fidelity, we set a normal prior with mean=0 and SD=2 to cover the range of most logistic regression coefficients [@Gelman2008]. 
  
  ```{r normalized degree interaction prior, eval=F, include=F}
  # As we did for the main effect for the scaled normalized degree of resources and consumers, we can use the same logic to generate a prior expectation for the statistical interaction between scaled normalized degrees.
  ## Resource species
  # probability of average species interacting
  pMean <- mean(full.df$r_ND) # arbitrarily chose r_ND, mean is the same as c_ND
  # probability of a relatively generalized resource interacting
  pGen_r <- pMean + sd(full.df$r_ND) 
  # probability of a relatively generalized consumer interacting
  pGen_c <- pMean + sd(full.df$c_ND) 
  # probability after statistical interaction
  pGen_rxc <- pGen_r * pGen_c
  # logistic regression coefficient, which is the difference in log-odds 
  # for a 1 SD unit increase in a species normalized degree.
  (mu_sc.r_NDxsc.c_ND <- log(pGen_rxc / (1-pGen_rxc)) - log(pMean / (1-pMean))) 
  ```
  
  
  ```{r simulate prior distribution for interaction between normalized, echo=FALSE, fig.cap="Simulations for four different choices of variance for a normal prior with mean=-0.06 for the statistical interaction between scaled normalized degree of resources and consumers.", eval=F, include=F}
  # This suggests that our prior expectation for the statistical interaction between scaled normalized degree of resources and consumers, if anything, will reduce the probability of an interaction. Note that this assumes we have already accounted for the main effects of each species scaled normalized degree.
  samp <- 1000
  sd_sc.r_NDxsc.c_ND <- data.frame(density = c(rnorm(n = samp, mean = mu_sc.r_NDxsc.c_ND, sd = 0.5),
                                               rnorm(n = samp, mean = mu_sc.r_NDxsc.c_ND, sd = 1),
                                               rnorm(n = samp, mean = mu_sc.r_NDxsc.c_ND, sd = 2),
                                               rnorm(n = samp, mean = mu_sc.r_NDxsc.c_ND, sd = 4)),
                                   standard_deviation = c(rep("sd_0.5", samp),
                                                          rep("sd_1", samp),
                                                          rep("sd_2", samp),
                                                          rep("sd_4", samp)))
  ggplot(sd_sc.r_NDxsc.c_ND, aes(x = density, fill = standard_deviation)) + 
    geom_density(alpha = 0.5) + 
    geom_vline(xintercept = mu_sc.r_NDxsc.c_ND, linetype = "dotted")
  # Since most standardized logistic regression coefficients are less than 5 [@Gelman2008], we feel that specifying a standard deviation of 2 is reasonable since it covers the range from 0 to $\pm$ 5. Therefore, we used the following prior: normal(mean = `r round(mu_sc.r_NDxsc.c_ND, 2)`, sd = 2).
  #Note that interpreting coefficients for statistical interactions between continuous variables is unwise without also plotting these effects. This is because, on the probability scale, the values of all variables in the model matter [@UCLA].
  ```
  
  
  ### Effect of mutualism on partner (a)symmetry in normalized degrees: **type_M.vs.A:sc.r_ND:sc.c_ND**
  
  We again set a normal prior with mean=0 and SD=2, because we had no strong prior expectation for how interaction type would modify the effect of (a)symmetry in normalized degree on partner fidelity. 
  
  ### Main effects and statistical interactions with observed interactions (sc.log.sum_)
  
  As with normalized degree, we expect the log of observed interactions for a resource or consumer to increase the probability of observing an interaction between co-occuring species. Therefore, we specified the same priors as we did for normalized degree for both main effects (normal(mean=1, SD=2) for **sc.log.sum_r, sc.log.sum_c**) and statistical interactions (normal(mean=0, SD=2) for **sc.log.sum_r:sc.log.sum_c, type_M.vs.A:sc.log.sum_r, type_M.vs.A:sc.log.sum_c, type_M.vs.A:sc.log.sum_r:sc.log.sum_c**).
  
  ### Random effects: resource_sp, consumer_sp, id_pair, network_id
  
  We specified a half-normal prior with mean=0 and SD=2 for each of the random effects in our model. This prior assumes that the variance in these effect sizes is small, but with sufficient evidence, the model can still estimate larger effects. Note that we specify this as a normal distribution in the code, but the R package we use constrains this prior to only positive values since standard deviations cannot be less than zero.
  
  \ 
  
  ## Model Analyses
  
  We used the *brms* package [@Burkner2017] and its default sampling behavior (four sampling chains for 2000 iterations each, discarding the first 1000 iterations as burn-in) to fit these statistical models. To mitigate bias in the posterior sampling distribution, we set **adapt_delta = 0.99** and **max_treedepth = 20**.
  
  For the majority of our inferences below, we keep coefficients on the log odds scale. This is because effect sizes are invariant on this scale. For intercept terms, we make inferences on the probability scale, and note the specific values at which this effect is measured. To calculate probabilities, we applied the inverse logit, $\frac{\exp(\beta)}{\exp(\beta)+1}$, to estimates reported in each plot. 
  
  ### Full model
  
  ```{r run model on full dataset, cache=TRUE, results="hide", message=F}
  full.brm <- brm(
    formula = interaction.formula, data = full.df, 
    prior = interaction.priors, algorithm = "sampling", 
    chains = 4, iter = 2000, warmup = 1000, 
    control = list(adapt_delta = 0.99, max_treedepth = 20))
  ```
  
  ```{r full-table, echo=FALSE, include=FALSE}
  tidy_full.brm <- tidy(full.brm, par_type = "non-varying", intervals = TRUE, prob = 0.95) 
  tidy_full.brm %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of fixed effects from our full model.", booktabs = T, linesep = "") %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r full-table-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of fixed effects from our full model.", out.width="75%", fig.align="center"}
  tidy_full.brm$term.ord <- factor(tidy_full.brm$term, 
                                   levels = c("type_M.vs.A:sc.r_ND:sc.c_ND", "sc.r_ND:sc.c_ND", "type_M.vs.A:sc.c_ND", "type_M.vs.A:sc.r_ND", "sc.c_ND", "sc.r_ND", "subtype_Poll.vs.Disp", "subtype_Herb.vs.Para", "type_M.vs.A", "Intercept"),
                                   #labels = c("Type:Resource(ND):Consumer(ND)","Resource(ND):Consumer(ND)","Type:Consumer(ND)","Type:Resource(ND)","Consumer(ND)","Resource(ND)","Subtype(Poll-vs-Disp)","Subtype(Herb-vs-Para)","Type(M-vs-A)","Intercept"), 
                                   ordered = T)
  full_table_plot <- tidy_full.brm %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey", size = 1.5) +
    geom_point(size = 3) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Coefficient estimate") +
    xlab("")
  full_table_plot
  ggsave("plots/full_model_plot.svg")
  ```
  
  \ 
  
  Our model provided a good fit to our data, explaining 40% of the variance.
  
  ```{r bayesian R2}
  bayes_R2(full.brm)
  ```
  
  \ 
  
  Below, we give a biological interpretation of each term in this model (fig. \ref{fig:full-table-plot}):
    
    - On the probability scale, the average level of partner fidelity across interaction types equals `r round(boot::inv.logit(tidy_full.brm$estimate[1]), 2)`($=\frac{\exp(\beta)}{\exp(\beta)+1}$, where $\beta=$ `r round(tidy_full.brm$estimate[1], 2)`) for the average level of normalized degree for a resource and consumer. 
  
  - Mutualistic interactions increase the log odds of partner fidelity by `r round(tidy_full.brm$estimate[2], 2)` relative to antagonistic interactions.
  
  - There is no clear difference between herbivory and parasitism on partner fidelity (95% credible intervals overlap with zero).
  
  - Pollination increases the log odds of partner fidelity by `r round(tidy_full.brm$estimate[4], 2)` relative to seed dispersers.
  
  - A 1 SD increase in resource normalized degree increases the log odds of partner fidelity by `r round(tidy_full.brm$estimate[5], 2)` across interaction types.
  
  - A 1 SD increase in consumer normalized degree increases the log odds of partner fidelity by `r round(tidy_full.brm$estimate[6], 2)` across interaction types.
  
  - Mutualistic interactions increase the positive effect of resource normalized degree on the log odds of partner fidelity by `r round(tidy_full.brm$estimate[7], 2)` relative to antagonistic interactions.
  
  - There is no clear evidence that mutualistic interactions (relative to antagonistic) modify the positive effect of consumer normalized degree on partner fidelity (95% credible intervals overlap with zero).
  
  - Across interaction types, symmetry between partners in their normalized degree increases the log odds of partner fidelity by `r round(tidy_full.brm$estimate[9], 2)`.
  
  - Mutualistic interactions (relative to antagonistic) decrease the effects of symmetry on the log odds of partner fidelity by `r round(tidy_full.brm$estimate[10], 2)`. Put another way, mutualistic interactions increase the effects of *asymmetry* on the log odds of partner fidelity by `r -1*round(tidy_full.brm$estimate[10], 2)`.
  
  
  ```{r explaining three-way interaction, include=FALSE}
  # these values are reported in the main text of the results 
  # to help explain the three-way statistical interaction
  leftpanel <- data.frame(
    type = c("M", "A"),
    subtype_Poll.vs.Disp = c(0, 0),
    subtype_Herb.vs.Para = c(0, 0),
    sc.c_ND = c(-1, -1),
    sc.r_ND = c(1, 1)
  )
  fitted(full.brm, newdata = leftpanel, scale = "response", re_formula = NA)
  rightpanel <- data.frame(
    type = c("M", "A"),
    subtype_Poll.vs.Disp = c(0, 0), 
    subtype_Herb.vs.Para = c(0, 0), 
    sc.c_ND = c(1, 1),
    sc.r_ND = c(-1, -1)
  )
  fitted(full.brm, newdata = rightpanel, scale = "response", re_formula = NA)
  ```
  
  
  ```{r three-way interaction plot, echo=FALSE, fig.cap="Asymmetric normalized degrees enhance partner fidelity in mutualistic vs. antagonistic interactions. Lines and bands correspond to predicted means and 95% credible intervals from our full model. Orange and blue colors correspond to mutualistic and antagonistic interactions, respectively. As expected, increasing normalized degree of resources enhances partner fidelity for both mutualistic and antagonistic interactions; however, mutualistic interactions have relatively higher partner fidelity when there is an asymmetry in partner normalized degrees. For example, mutualistic interactions enhance partner fidelity for specialized consumers (scaled ND = -1) and generalist resources (scaled ND > 0) relative to antagonistic interactions. Likewise, mutualistic interactions enhance partner fidelity for generalized consumers (scaled ND = 1) and specialist resources (scaled ND < 0) relative to antagonistic interactions.", include=FALSE}
  # we like this perspective the most.
  # as above, suggests to me that specialized consumers have a much higher probability of interacting with generalist resources 
  # in mutualistic vs. antagonistic interactions
  # color-blind friendly palette from: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  # cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # not good for gray-scale printouts
  AvM_colors <- c("#d95f02","#1b9e77") # chose new colors using ColorBrewer http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3 that are colorblind and print friendly, for a qualitative data display # cbPalette[6:7]
  asymmetry_SvG <- plot(marginal_effects(full.brm, effects = "sc.r_ND:type",
                                         conditions = data.frame(sc.c_ND = c(-1,1),
                                                                 row.names = c("Specialist consumer\n(scaled ND = -1)", "Generalist consumer\n(scaled ND = 1)")),
                                         probs = c(0.025,0.975)))[[1]] +
    ggplot2::ylab("Partner fidelity") +
    ggplot2::xlab("Resource Normalized Degree (scaled)") +
    ggplot2::scale_color_manual(values = AvM_colors, name = "Interaction type", labels = c("Antagonistic", "Mutualistic")) +
    ggplot2::scale_fill_manual(values = AvM_colors, name = "Interaction type", labels = c("Antagonistic", "Mutualistic")) +
    # not working # ggplot2::scale_linetype_manual(values = c("solid","dashed"), name = "Interaction type", labels = c("Antagonistic", "Mutualistic")) +
    ggplot2::theme(strip.background = element_blank())
  asymmetry_SvG
  ggsave("plots/asymmetry_specialist_vs_generalist.svg")
  ```
  
  In addition to these fixed effects, there is considerable variation in partner fidelity explained by our random effects. In particular, **network_id** has a strong effect compared to the other sources of variation (fig. \ref{fig:full-table-random-plot}).
  
  ```{r full-table-random, echo=FALSE, include=FALSE}
  tidy(full.brm, par_type = "hierarchical", intervals = TRUE, prob = 0.95) %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of random effects from our full model.", booktabs = T) %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r full-table-random-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of random effects from our full model.", out.width="75%", fig.align="center"}
  tidy_full.brm.random <- tidy(full.brm, par_type = "hierarchical", intervals = TRUE, prob = 0.95)
  tidy_full.brm.random$term.ord <- factor(tidy_full.brm.random$term, 
                                          levels = c("sd_consumer_sp__Intercept", "sd_resource_sp__Intercept", "sd_id_pair__Intercept", "sd_network_id__Intercept"),
                                          labels = c("Consumer ID","Resource ID","Pair ID","Network ID"), ordered = T)
  tidy_full.brm.random %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey", size = 1.5) +
    geom_point(size = 3) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Standard Deviations") +
    xlab("")
  ```
  
  \ 
  
  While the model above is useful for determining the effects of mutualistic relative to antagonistic interactions on partner fidelity, it is also useful to separately estimate relationships for mutualistic and antagonistic interactions. We can do this by removing the intercept term from the model, and estimating unique relationships for normalized degree for both mutualistic and antagonistic interactions.
  
  ```{r run no intercept model on full dataset, cache=TRUE, results="hide", message=F}
  no.intercept_formula <- brmsformula(
    connected ~ -1 + type + subtype_Herb.vs.Para + subtype_Poll.vs.Disp + 
      type:sc.r_ND + type:sc.c_ND + 
      type:sc.r_ND:sc.c_ND +
      (1 | resource_sp) + (1 | consumer_sp) + (1 | id_pair) + (1 | network_id),
    family = bernoulli(link = "logit"))
  no.intercept_priors <- c(
    set_prior("normal(0,2)", class = "b", coef = "typeA"), # corresponds to probability of 0.5
    set_prior("normal(0,2)", class = "b", coef = "typeM"), # corresponds to probability of 0.5
    set_prior("normal(0,2)", class = "b", coef = "subtype_Herb.vs.Para"),
    set_prior("normal(0,2)", class = "b", coef = "subtype_Poll.vs.Disp"),
    set_prior("normal(1,2)", class = "b", coef = "typeA:sc.r_ND"),
    set_prior("normal(1,2)", class = "b", coef = "typeM:sc.r_ND"),
    set_prior("normal(1,2)", class = "b", coef = "typeA:sc.c_ND"),
    set_prior("normal(1,2)", class = "b", coef = "typeM:sc.c_ND"),
    set_prior("normal(0,2)", class = "b", coef = "typeA:sc.r_ND:sc.c_ND"),
    set_prior("normal(0,2)", class = "b", coef = "typeM:sc.r_ND:sc.c_ND"),
    set_prior("normal(0,2)", class = "sd"))
  no.intercept_brm <- brm(
    formula = no.intercept_formula, data = full.df, 
    prior = no.intercept_priors, algorithm = "sampling", 
    chains = 4, iter = 2000, warmup = 1000, 
    control = list(adapt_delta = 0.99, max_treedepth = 20))
  ```
  
  \ 
  
  ```{r no-intercept-table, echo=FALSE, include=FALSE}
  tidy_no.intercept_brm <- tidy(no.intercept_brm, par_type = "non-varying", intervals = TRUE, prob = 0.95) 
  tidy_no.intercept_brm %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of fixed effects from model without intercepts.", booktabs = T, linesep = "") %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r no-intercept-table-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of fixed effects from model without intercepts.", out.width="75%", fig.align="center"}
  tidy_no.intercept_brm$term.ord <- factor(tidy_no.intercept_brm$term, 
                                           levels = c("typeM:sc.r_ND:sc.c_ND", "typeA:sc.r_ND:sc.c_ND", "typeM:sc.c_ND", "typeA:sc.c_ND", "typeM:sc.r_ND", "typeA:sc.r_ND", "subtype_Poll.vs.Disp", "subtype_Herb.vs.Para", "typeM", "typeA"),
                                           #labels = c("Type:Resource(ND):Consumer(ND)","Resource(ND):Consumer(ND)","Type:Consumer(ND)","Type:Resource(ND)","Consumer(ND)","Resource(ND)","Subtype(Poll-vs-Disp)","Subtype(Herb-vs-Para)","Type(M-vs-A)","Intercept"), 
                                           ordered = T)
  tidy_no.intercept_brm %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey", size = 1.5) +
    geom_point(size = 3) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Coefficient estimate") +
    xlab("")
  ```
  
  \ 
  
  The results from this model (fig. \ref{fig:no-intercept-table-plot}) can be interpreted as follows:
    
    - The average level of partner fidelity for antagonistic interactions is `r round(boot::inv.logit(tidy_no.intercept_brm$estimate[1]), 2)` ($=\frac{\exp(\beta)}{\exp(\beta)+1}$, where $\beta=$ `r round(tidy_no.intercept_brm$estimate[1], 2)`) at the average level of normalized degree for a resource and consumer.
  
  - The average probability of partner fidelity for mutualistic interactions is `r round(boot::inv.logit(tidy_no.intercept_brm$estimate[2]), 2)` ($=\frac{\exp(\beta)}{\exp(\beta)+1}$, where $\beta=$ `r round(tidy_no.intercept_brm$estimate[2], 2)`)  at the average level of normalized degree for a resource and consumer.
  
  - A 1 SD increase in resource normalized degree in antagonistic interactions increases the log odds of an interaction by `r round(tidy_no.intercept_brm$estimate[5], 2)`.
  
  - A 1 SD increase in resource normalized degree in mutualistic interactions increases the log odds of an interaction by `r round(tidy_no.intercept_brm$estimate[6], 2)`.
  
  - A 1 SD increase in consumer normalized degree in antagonistic interactions increases the log odds of an interaction by `r round(tidy_no.intercept_brm$estimate[7], 2)`.
  
  - A 1 SD increase in consumer normalized degree in mutualistic interactions increases the log odds of an interaction by `r round(tidy_no.intercept_brm$estimate[8], 2)`.
  
  - Symmetry in partner normalized degrees increases the log odds of partner fidelity by `r round(tidy_no.intercept_brm$estimate[9], 2)` for antagonistic interactions.
  
  - Asymmetry in partner normalized degrees increases the log odds of partner fidelity by `r -1*round(tidy_no.intercept_brm$estimate[10], 2)` for mutualistic interactions.
  
  
  ### Analysis of data subset
  
  Before testing for potential biases due to sampling effort, we test whether we observe the same results for the data subset.
  
  ```{r run weighted model, cache=TRUE, results="hide", message=F}
  weighted_subset.brm <- brm(
    formula = interaction_subset.formula, data = weighted_subset.df, 
    prior = interaction_subset.priors, algorithm = "sampling", 
    chains = 4, iter = 2000, warmup = 1000,
    control = list(adapt_delta = 0.99, max_treedepth = 20))
  ```
  
  ```{r weighted-table, echo=FALSE, include=FALSE}
  tidy(weighted_subset.brm, par_type = "non-varying", intervals = TRUE, prob = 0.95) %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of fixed effects from our weighted subset model.", booktabs = T) %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r weighted-table-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of fixed effects from our weighted subset model.", out.width="75%", fig.align="center"}
  tidy_weighted_subset_brm <- tidy(weighted_subset.brm, par_type = "non-varying", intervals = TRUE, prob = 0.95)
  tidy_weighted_subset_brm$term.ord <- factor(tidy_weighted_subset_brm$term, 
                                              levels = c("type_M.vs.A:sc.r_ND:sc.c_ND", "sc.r_ND:sc.c_ND", "type_M.vs.A:sc.c_ND", "type_M.vs.A:sc.r_ND", "sc.c_ND",
                                                         "sc.r_ND", "subtype_Poll.vs.Disp", # not present # "subtype_Herb.vs.Para", 
                                                         "type_M.vs.A", "Intercept"),
                                              #labels = c("Type:Resource(ND):Consumer(ND)","Resource(ND):Consumer(ND)","Type:Consumer(ND)","Type:Resource(ND)","Consumer(ND)","Resource(ND)","Subtype(Poll-vs-Disp)","Subtype(Herb-vs-Para)","Type(M-vs-A)","Intercept"), 
                                              ordered = T)
  tidy_weighted_subset_brm %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey", size = 1.5) +
    geom_point(size = 3) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Coefficient estimate") +
    xlab("")
  ```
  
  Note that the **Intercept** and effect of **type_M.vs.A** is higher in this subset of data (fig. \ref{fig:weighted-table-plot}). This is likely due to the fact that herbivory had a tendency to increase partner fidelity (**subtype_Herb.vs.Para** in fig. \ref{fig:full-table-plot}); therefore, its absence from the data subset enhances the average estimate of partner fidelity as well as the apparent effect of mutualistic interactions. Other terms are similar, except that there is less clear evidence of mutualistic interactions enhancing the effect of resource normalized degree (**type_M.vs.A:sc.r_ND** in fig. \ref{fig:weighted-table-plot}), and less clear evidence of symmetry in partner normalized degrees to positively affect partner fidelity (**sc.r_ND:sc.c_ND** in fig. \ref{fig:weighted-table-plot}). Importantly, there is still clear evidence that mutualistic interactions generally enhance partner fidelity (**type_M.vs.A** in fig. \ref{fig:weighted-table-plot}) and also enhance the effects of *asymmetry* in partner normalized degrees on partner fidelity (**type_M.vs.A:sc.r_ND:sc.c_ND** in fig. \ref{fig:weighted-table-plot}).
  
  As with the full model, we observe a particularly strong effect of **network_id** compared to the other sources of variation (fig. \ref{fig:weighted-table-random-plot}).
  
  ```{r weighted-table-random, echo=FALSE, include=FALSE}
  tidy(weighted_subset.brm, par_type = "hierarchical", intervals = TRUE, prob = 0.95) %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of random effects from our weighted subset model.", booktabs = T) %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r weighted-table-random-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of random effects from our weighted subset model.", out.width = "75%", fig.align="center"}
  tidy_weighted_subset_brm.random <- tidy(weighted_subset.brm, par_type = "hierarchical", intervals = TRUE, prob = 0.95)
  tidy_weighted_subset_brm.random$term.ord <- factor(tidy_weighted_subset_brm.random$term, 
                                                     levels = c("sd_consumer_sp__Intercept", "sd_resource_sp__Intercept", "sd_id_pair__Intercept",
                                                                "sd_network_id__Intercept"),
                                                     labels = c("Consumer ID","Resource ID","Pair ID","Network ID"), 
                                                     ordered = T)
  tidy_full.brm.random %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey", size = 1.5) +
    geom_point(size = 3) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Standard Deviations") +
    xlab("")
  ```
  
  ### Testing for bias due to sampling effort
  
  Our goal here is to test whether variation in sampling effort between mutualistic and antagonistic networks could explain the results we observed in the previous models, in particular, the three-way statistical interaction with resource and consumer normalized degree. 
  
  Using the subset of weighted data, we first create a long dataset where each row corresponds to an observed interaction.
  
  ```{r long dataset, eval=FALSE}
  # focus dataset
  weighted_long <- weighted_subset.df %>%
    select(network_id, type, subtype, resource_sp, consumer_sp, connected, sum_shared)
  # for each interaction, replicate it as many times as it was observed
  weighted_long.list <- list()
  for(i in 1:nrow(weighted_long)){ 
    if(weighted_long$sum_shared[i] > 0){
      weighted_long.list[[i]] <- slice(weighted_long[i, ], 
                                       rep(1:n(), 
                                           each = weighted_long$sum_shared[i]))
    } else {
      weighted_long.list[[i]] <- weighted_long[i, ]
    }
  }
  # tidy into a dataframe
  weighted_long.df <- plyr::ldply(weighted_long.list) %>% 
    select(-sum_shared)
  ```
  
  \ 
  
  To control for sampling effort, we conducted the following procedure:
    
    1. Pair each mutualistic network with a random subset of the antagonistic networks;
  2. Randomly subsample each network down to the interaction count of the paired network with the least number of interactions sampled;
  3. Estimate the three-way interaction term for the subsampled data;
  4. Repeat this procedure 1,000 times;
  5. Compare the distribution of estimates to the reference estimated from the subset of networks with weighted data (i.e. point representing **type_M.vs.A:sc.r_ND:sc.c_ND** in fig. \ref{fig:weighted-table-plot}).
  
  ```{r subsampling code, eval=FALSE}
  n_sims <- 1000
  three_way.FE <- c() # vector of fixed-effect estimate for three-way interaction term
  three_way.P <- c() # vector of P-values for the three-way term
  for(i in 1:n_sims){
    # subsample Antagonistic networks to match the number of mutualistic ones
    subsample_A_networks <- distinct(filter(weighted_long.df, type == "A"), network_id) %>% 
      sample_n(size = 19, replace = F) 
    
    # no need to subsample because there are fewer
    M_networks <- distinct(filter(weighted_long.df, type == "M"), network_id) 
    
    # randomly pair them together (note randomization was already done with subsample_A_networks)
    random_pairs <- data.frame(network_pairs = LETTERS[1:19],
                               M = as.character(M_networks$network_id),
                               A = as.character(subsample_A_networks$network_id))
    
    # create data frame to determine subsample size for each network
    subsample_df <- left_join(random_pairs, 
                              select(weighted_subset.df.network_level, 
                                     M = network_id, M_sum = sum_shared)) %>%
      left_join(., select(weighted_subset.df.network_level, 
                          A = network_id, A_sum = sum_shared)) %>%
      # subsample is based on smallest number of observed interactions in a network
      mutate(subsample_size = ifelse(M_sum < A_sum, M_sum, A_sum)) %>%
      select(network_pairs, M, A, subsample_size) %>%
      # organize for joining to other datasets
      gather(type, network_id, -subsample_size, -network_pairs)
    
    # now I need to make the new data frame based on the random subsample
    new_data <- weighted_long.df %>%
      # filter to randomly sampled networks
      filter(network_id %in% c(as.character(M_networks$network_id), 
                               as.character(subsample_A_networks$network_id))) %>%
      # following example from: 
      # https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
      # to get different sample sizes for each group
      group_by(network_id) %>%
      nest() %>%
      left_join(., subsample_df) %>%
      mutate(samp = map2(data, subsample_size, sample_n)) %>%
      select(-data) %>%
      unnest(samp) %>%
      # now that I've subsampled, I can reduce the dataset to the normal size for analysis
      # i.e. only 1 interaction per pair per network
      distinct %>%
      # add back in some data for the analysis
      # note that sc.r_ND and sc.c_ND are calculated from the original data and 
      # not the subsampled interactions. this is because their normalized degree 
      # is based on the entire interaction network rather than their partner fidelity 
      # across networks
      left_join(., select(weighted_subset.df, 
                          network_id, resource_sp, consumer_sp, id_pair, sc.r_ND, sc.c_ND))
    
    # run test with glmer
    # Note that we use a frequentist, rather than Bayesian approach, to calculate 
    # the three-way interaction term in order to speed up the simulation, but this should not 
    # matter since our results are not affected by our choice of priors 
    # (see **Robustness to Different Priors** section).
    new_glmer <- glmer(connected ~ subtype + type*sc.r_ND*sc.c_ND + 
                         (1|consumer_sp) + (1|resource_sp) + (1|id_pair) + (1|network_id), 
                       data = new_data, family = binomial(link = "logit"), 
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    
    # get model effects for 3-way interaction term
    new_glmer_effects <- broom::tidy(new_glmer, effects = "fixed")
    three_way.FE[i] <- filter(new_glmer_effects, term == "typeM:sc.r_ND:sc.c_ND")$estimate
    three_way.P[i] <- filter(new_glmer_effects, term == "typeM:sc.r_ND:sc.c_ND")$p.value
    # print(paste(i/n_sims*100,"% done!")) # track progress
  }
  # write data from simulation 
  write_csv(x = data.frame(FixedEffect = three_way.FE, P = three_way.P), 
            path = "data/three_way_permutation_data.csv")
  ```
  
  
  ```{r plot subsample of three-way term, echo=FALSE, results="hide", message=F}
  three_way_perms <- read_csv("data/three_way_permutation_data.csv")
  # reference statistic based on three-way term in Figure 7 of supplementary material
  ref.FE <- filter(tidy_weighted_subset_brm, term == "type_M.vs.A:sc.r_ND:sc.c_ND")$estimate #-0.68
  # looks to me like equalizing the sample size doesn't alter the distribution of the test statistic,
  # suggesting that our analysis is robust
  ggplot(three_way_perms, aes(x = FixedEffect)) +
    geom_histogram(alpha = 0.5, color = "grey") +
    geom_vline(xintercept = ref.FE, linetype = "dotted") +
    xlab("Coefficient estimate") +
    ylab("Count")
  # proportion of times that equalizing sample size reduces effect size
  sum(three_way_perms$FixedEffect < ref.FE) / length(three_way_perms$FixedEffect) # 55% of the time.
  # i.e. since the reference statistic is not larger than 90% or more of the subsampled
  # distribution values, then sampling intensity does not matter for our results.
  # proportion of times that equalizing sample size flips the sign of the effect
  sum(three_way_perms$FixedEffect > 0) / length(three_way_perms$FixedEffect) # 2%
  # proportion of times the 3-way interaction is statistically significant
  sum(three_way_perms$P < 0.05) / length(three_way_perms$P) # 35%
  ```
  
  \ 
  
  If variation in sampling effort between mutualistic and antagonistic networks were driving our results, then we should expect our reference estimate for the three-way interaction term to be lower than 90% or more of the subsampled values (i.e. a conservative one-tailed permutation test). In contrast, we see that the reference is lower than only `r round(sum(three_way_perms$FixedEffect > ref.FE) / length(three_way_perms$FixedEffect),2)*100`% of the subsamples. This indicates that variation in sampling effort between mutualistic and antagonistic networks cannot explain our results. 
  
  ## Robustness to Different Priors
  
  To test how robust our results were to a different choice of priors, we re-analyzed the full model with the default priors from the *brms* package, which are designed to be non or very weakly informative. Note that we had to increase **adapt_delta = 0.999** to avoid bias in the posterior sampling of this model.
  
  ```{r run no priors model on full dataset, cache=TRUE, results="hide", message=F}
  nopriors.brm <- brm(
    formula = interaction.formula, data = full.df, 
    # prior = interaction.priors, # remove and use default priors
    algorithm = "sampling", chains = 4, iter = 2000, warmup = 1000, 
    control = list(adapt_delta = 0.999, max_treedepth = 20))
  ```
  
  
  ```{r nopriors-table, echo=FALSE, include=FALSE}
  tidy(nopriors.brm, par_type = "non-varying", intervals = TRUE, prob = 0.95) %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of fixed effects from the no priors model.", booktabs = T) %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r nopriors-table-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of fixed effects from our full model.", out.width="75%", fig.align="center"}
  tidy_nopriors.brm <- tidy(nopriors.brm, par_type = "non-varying", intervals = TRUE, prob = 0.95)
  tidy_nopriors.brm$term.ord <- factor(tidy_nopriors.brm$term, 
                                       levels = c("type_M.vs.A:sc.r_ND:sc.c_ND", "sc.r_ND:sc.c_ND", "type_M.vs.A:sc.c_ND", "type_M.vs.A:sc.r_ND", "sc.c_ND", "sc.r_ND", "subtype_Poll.vs.Disp", "subtype_Herb.vs.Para", "type_M.vs.A", "Intercept"),
                                       #labels = c("Type:Resource(ND):Consumer(ND)","Resource(ND):Consumer(ND)","Type:Consumer(ND)","Type:Resource(ND)","Consumer(ND)","Resource(ND)","Subtype(Poll-vs-Disp)","Subtype(Herb-vs-Para)","Type(M-vs-A)","Intercept"), 
                                       ordered = T)
  tidy_nopriors.brm %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey") +
    geom_point() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Coefficient estimate") +
    xlab("")
  ```
  
  ```{r nopriors-table-random, echo=FALSE, include=FALSE}
  tidy(nopriors.brm, par_type = "hierarchical", intervals = TRUE, prob = 0.95) %>%
    transmute(Term = term,
              Estimate = round(estimate, 2),
              `2.5%` = round(lower, 2),
              `97.5%` = round(upper, 2)) %>%
    knitr::kable(., caption = "Mean and 95\\% credible intervals of random effects from the no priors model.", booktabs = T) %>%
    kable_styling(latex_options = c("striped", "hold_position"))
  ```
  
  ```{r nopriors-table-random-plot, echo=FALSE, fig.cap = "Mean and 95\\% credible intervals of random effects from the no priors model.", out.width = "75%", fig.align="center"}
  tidy_nopriors_brm.random <- tidy(nopriors.brm, par_type = "hierarchical", intervals = TRUE, prob = 0.95)
  tidy_nopriors_brm.random$term.ord <- factor(tidy_nopriors_brm.random$term, 
                                              levels = c("sd_consumer_sp__Intercept", "sd_resource_sp__Intercept", "sd_id_pair__Intercept",
                                                         "sd_network_id__Intercept"),
                                              labels = c("Consumer ID","Resource ID","Pair ID","Network ID"), 
                                              ordered = T)
  tidy_nopriors_brm.random %>%
    ggplot(., aes(x = term.ord, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "grey") +
    geom_point() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Standard Deviations") +
    xlab("")
  ```
  
  \
  