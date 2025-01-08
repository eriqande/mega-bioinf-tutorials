Using List Columns (tibble, tidyr, purrr)
================

## Introduction

This is a tutorial intended for my esteemed labmates and members of the
weekly MEGA-bioinformatics group. I have been intending to do a short
tutorial on this way of wrangling data, and now we have the perfect use
case for exploring it.

Our goal today is to use `tidyr::nest()` and functions from the `purrr`
package to streamline some analyses that we wish to do on various
subsets of the data. Our use case is Laura’s linear-mixed modeling on
steelhead migration timing in the Central Valley as a function of
certain GREB1L loci.

First, let’s download Laura’s notebook and investigate it to see what we
are up against. You can go to our google drive folder (you have access
with your email that the Mega-bioinf calendar notifications come to) at
<https://drive.google.com/drive/folders/1M0VeY-kwtuWu3lJjX8COOEoNeF5mrzXy?usp=sharing>.

Open up that folder and download `Associationtests.html`, then open it
from your downloads folder.

The action that we want to be dealing with is in the **Running Models**
section of the notebook. What Laura is doing is running similar analyses
on different subsets (different loci and different populations) of a big
data set. These analyses take a data frame as input, and they produce
complex output, so they are not amenable to a quick `summarise()` on a
grouped data frame. But they are amenable to operation on nested data
frames. The basic steps in the analysis are:

1.  Make a tibble that has a subset of the data (subset by locus and
    population)
2.  Make a design matrix that is specific to the locus. The main inputs
    that change for these are the genotypes.
3.  Add the design matrix columns on by joining on genotype
4.  Run a series of linear mixed models on each subset using the
    `lmer()` function. The different models are:
    - `mo_da ~ I(d_add) + I(d_dom_with_s) + (1|year)`
    - `mo_da ~ I(d_add) + I(d_dom_with_s) + sex + (1|year)`
    - `mo_da ~ I(d_add) + I(d_dom_with_s) + age_spawn + (1|year)`
5.  We would like to be able to access a summary of each of those
    models, and maybe even make some plots, programmatically of the
    results.

Laura put this all together nicely, and ended up running each of the
models on each of the subsets by hand. This is great for exploratory
work; however, when you find that you have copied and pasted the same
block of code, and then changed a few variable names in it more than two
or three times, it is time to start thinking about abstracting those
steps into a function that can be applied to each different subset of
the data. This leads to less code to maintain—if you want to make a
change in your code you don’t have to make it in 15 different places—and
also lessens the chance of making a typo or other error. That said,
sometimes it is hard to follow things back through a chain of functions,
but, for the most part, it is good to break repetitive things into
functions, and it also works well with purr.

### Packages

We need to have a few packages here. The following code will install
necessary packages if you don’t already have them.

``` r
need_em <- c("tidyverse", "lme4", "lubridate", "lmerTest", "cowplot")
please_install <- setdiff(need_em, rownames(installed.packages()))
if(length(please_install) > 0) {
  install.packages(please_install)
}
```

Once those are installed, we can load them:

``` r
library(tidyverse)
library(lme4)
library(lmerTest)
library(lubridate)
library(cowplot)
```

### Getting our data and some early processing

From the same google drive folder linked above, you can download the
data set, `fyke_grebs.rds` and put it in the current working directory.

Then read it in and have a look at it:

``` r
fyke_grebs <- read_rds("fyke_grebs.rds")
```

## Make a tibble of tibbles—each with a subet of data

With this step, we are going to
