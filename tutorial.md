Using List Columns (tibble, tidyr, purrr)
================

- [Introduction](#introduction)
  - [Packages](#packages)
  - [Getting our data and some early
    processing](#getting-our-data-and-some-early-processing)
- [Make a tibble of tibbles—each with a subet of
  data](#make-a-tibble-of-tibbleseach-with-a-subet-of-data)
- [Use `map()` get get the genotypes at each
  locus](#use-map-get-get-the-genotypes-at-each-locus)
- [*Typed* forms of `map()`](#typed-forms-of-map)
  - [Number of genotypes using
    `map_int()`](#number-of-genotypes-using-map_int)
- [Dealing with three different
  models](#dealing-with-three-different-models)
- [A function to run `lmer()`](#a-function-to-run-lmer)
- [Using `pmap()` to do `run_lmer()` over all rows of the
  tibble](#using-pmap-to-do-run_lmer-over-all-rows-of-the-tibble)
- [One thing you might to see is the summary of all the
  models](#one-thing-you-might-to-see-is-the-summary-of-all-the-models)
- [Or use broom.mixed](#or-use-broommixed)
- [Let’s do a quick look at pvalues over loci and
  groups](#lets-do-a-quick-look-at-pvalues-over-loci-and-groups)

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
    also Nimbus Hatchery versus not Nimbus Hatchery)
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

fyke_grebs
```

    ## # A tibble: 17,520 × 24
    ##    indiv   spawner_group sex    year hatchery ArchiveID LENGTH Method collection
    ##    <chr>   <chr>         <chr> <dbl> <chr>    <chr>      <dbl> <chr>  <chr>     
    ##  1 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  2 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  3 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  4 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  5 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  6 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  7 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  8 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  9 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ## 10 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ## # ℹ 17,510 more rows
    ## # ℹ 15 more variables: SpawnYear <dbl>, kid_year <dbl>, kid_hatchery <chr>,
    ## #   ma_hatchery <chr>, age_spawn <dbl>, mo_dat <dbl>, mo_da <dbl>,
    ## #   program_assn <chr>, Locus <chr>, Locus_new <chr>, genecopy_1 <chr>,
    ## #   genecopy_2 <chr>, genotype_number <dbl>, genotype <chr>, typed_final <chr>

There is a lot of stuff there. Let’s look at the number of records for
different subsets we will be breaking it out into:

``` r
fyke_grebs %>%
  count(collection, Locus_new)
```

    ## # A tibble: 100 × 3
    ##    collection Locus_new                         n
    ##    <chr>      <chr>                         <int>
    ##  1 CH         Omy_Ch28_11622532_B_position1   472
    ##  2 CH         Omy_Ch28_11623412_B_position1   472
    ##  3 CH         Omy_Ch28_11625241_C_position1   469
    ##  4 CH         Omy_Ch28_11629294_D_position1   461
    ##  5 CH         Omy_Ch28_11641623_D_position1   470
    ##  6 CH         Omy_Ch28_11658853_F_position1   472
    ##  7 CH         Omy_Ch28_11667578_F_position1   471
    ##  8 CH         Omy_Ch28_11676622_G_position1   375
    ##  9 CH         Omy_Ch28_11678603_G_position1   472
    ## 10 CH         Omy_Ch28_11683310_H_position1   434
    ## # ℹ 90 more rows

We also count up states of variables used in the modeling.

Year:

``` r
fyke_grebs %>%
  count(year)
```

    ## # A tibble: 5 × 2
    ##    year     n
    ##   <dbl> <int>
    ## 1  2017   733
    ## 2  2018  3801
    ## 3  2019  5147
    ## 4  2020  4089
    ## 5  2021  3750

Sex:

``` r
fyke_grebs %>%
  count(sex)
```

    ## # A tibble: 3 × 2
    ##   sex        n
    ##   <chr>  <int>
    ## 1 ?        766
    ## 2 Female  9053
    ## 3 Male    7701

age_spawn:

``` r
fyke_grebs %>%
  count(age_spawn)
```

    ## # A tibble: 5 × 2
    ##   age_spawn     n
    ##       <dbl> <int>
    ## 1         1    69
    ## 2         2 10140
    ## 3         3  1999
    ## 4         4    24
    ## 5        NA  5288

OK, when we use these in modeling, we will want them to be factors.
Let’s make factor versions of each, noting that `?` is missing data for
sex.

And we also want to have another column called `group` that tells us
whether the ancestry of the fish is Central Valley (all collections
except NH) or not (NH)

``` r
fyke_grebs2 <- fyke_grebs %>%
  mutate(
    sex_f = ifelse(sex == "?", NA, sex) %>% factor(),
    year_f = factor(year),
    age_spawn_f = factor(age_spawn),
    group = ifelse(collection == "NH", "NH", "CV")
  )

# then check that this gave reasonable results
levels(fyke_grebs2$year_f)
```

    ## [1] "2017" "2018" "2019" "2020" "2021"

``` r
levels(fyke_grebs2$age_spawn_f)
```

    ## [1] "1" "2" "3" "4"

``` r
levels(fyke_grebs2$sex_f)
```

    ## [1] "Female" "Male"

OK! `fyke_grebs2` is, now, what we want to use.

## Make a tibble of tibbles—each with a subet of data

With this step, we are going to transform fyke_grebs2 into a new tibble
that has a list column that itself holds different tibbles. To do this
we use the `nest()` function from the ‘tidyr’ package, which is loaded
as part of the tidyverse. With no other arguments, `nest()` operates
like a `summarise()` on a *grouped tibble*. But the summary option is
simply *take all the rows within each level of grouping variables and
squash them down into a separate tibble in a list column that is named
`data` by default*. Let’s do it:

``` r
fg_nests <- fyke_grebs2 %>%
  group_by(group, Locus_new) %>%
  nest()

fg_nests
```

    ## # A tibble: 50 × 3
    ## # Groups:   group, Locus_new [50]
    ##    Locus_new                      group data               
    ##    <chr>                          <chr> <list>             
    ##  1 Omy_Ch28_greb1_mhap1_position1 CV    <tibble [631 × 26]>
    ##  2 Omy_Ch28_11658853_F_position1  CV    <tibble [637 × 26]>
    ##  3 Omy_Ch28_11622532_B_position1  CV    <tibble [638 × 26]>
    ##  4 Omy_Ch28_11623412_B_position1  CV    <tibble [637 × 26]>
    ##  5 Omy_Ch28_11641623_D_position1  CV    <tibble [636 × 26]>
    ##  6 Omy_Ch28_11629294_D_position1  CV    <tibble [626 × 26]>
    ##  7 Omy_Ch28_greb1_mhap2_position1 CV    <tibble [627 × 26]>
    ##  8 Omy_Ch28_11625241_C_position1  CV    <tibble [635 × 26]>
    ##  9 Omy_Ch28_greb1_mhap6_position1 CV    <tibble [637 × 26]>
    ## 10 Omy_Ch28_greb1_mhap8_position1 CV    <tibble [635 × 26]>
    ## # ℹ 40 more rows

Cool! We see that each row of this new tibble has a different
combination of `collection` and `Locus_new` and the default printing of
the list column `data` shows the size of each tibble that is an element
of the list.

This is all well and good, but how do we now get at portions of that
tibble (if we need them) in a tidyverse-like way? A simple `mutate()`
will not work, because `mutate()` expects a vectorized function, and not
many functions are vectorized to work over the elements of a list. Aha!
But, like `lapply()` the `map()` family of functions works over the
elements of a list, and we can use that inside of a `mutate()`. Next
section shows how.

## Use `map()` get get the genotypes at each locus

If we look back at Laura’s original notebook we see that one of the
inputs that we are going to need for each model run is the actual
genotypes at each locus. If we just do a find in the notebook for
`genotype = c(` we can see that the genotypes are always specified in a
vector of genotypes specified like `A/G` and always in sorted order.

Let’s make a new nested tibble where we group just on `Locus_new` so
that we are sure we have all the alleles seen in the whole data set.
Note that there are some rows of Locus_new that are NA, so we will toss
those.

``` r
fg_locus_nest <- fyke_grebs2 %>%
  filter(!is.na(Locus_new)) %>%
  group_by(Locus_new) %>%
  nest()

fg_locus_nest
```

    ## # A tibble: 24 × 2
    ## # Groups:   Locus_new [24]
    ##    Locus_new                      data               
    ##    <chr>                          <list>             
    ##  1 Omy_Ch28_greb1_mhap1_position1 <tibble [739 × 27]>
    ##  2 Omy_Ch28_11658853_F_position1  <tibble [745 × 27]>
    ##  3 Omy_Ch28_11622532_B_position1  <tibble [746 × 27]>
    ##  4 Omy_Ch28_11623412_B_position1  <tibble [745 × 27]>
    ##  5 Omy_Ch28_11641623_D_position1  <tibble [744 × 27]>
    ##  6 Omy_Ch28_11629294_D_position1  <tibble [732 × 27]>
    ##  7 Omy_Ch28_greb1_mhap2_position1 <tibble [735 × 27]>
    ##  8 Omy_Ch28_11625241_C_position1  <tibble [741 × 27]>
    ##  9 Omy_Ch28_greb1_mhap6_position1 <tibble [745 × 27]>
    ## 10 Omy_Ch28_greb1_mhap8_position1 <tibble [739 × 27]>
    ## # ℹ 14 more rows

We see that is 24 different SNPs. Let’s just look at the tibble in the
first row of the `data` column:

``` r
fg_locus_nest$data[[1]]
```

    ## # A tibble: 739 × 27
    ##    indiv   spawner_group sex    year hatchery ArchiveID LENGTH Method collection
    ##    <chr>   <chr>         <chr> <dbl> <chr>    <chr>      <dbl> <chr>  <chr>     
    ##  1 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  2 M096526 9/1/16        Fema…  2017 fyke     M096526      448 FYKE … FRH       
    ##  3 M096527 9/7/16        Fema…  2017 fyke     M096527      452 FYKE … CH        
    ##  4 M096529 9/11/16       Male   2017 fyke     M096529      445 FYKE … MRH       
    ##  5 M096530 9/13/16       Male   2017 fyke     M096530      380 FYKE … CH        
    ##  6 M096532 9/22/16       Fema…  2017 fyke     M096532      431 FYKE … CH        
    ##  7 M096533 9/25/16       Male   2017 fyke     M096533      480 FYKE … CH        
    ##  8 M096534 9/26/16       Male   2017 fyke     M096534      478 FYKE … CH        
    ##  9 M096535 10/1/16       Male   2017 fyke     M096535      500 FYKE … CH        
    ## 10 M096536 11/30/16      Male   2017 fyke     M096536      500 FYKE … CH        
    ## # ℹ 729 more rows
    ## # ℹ 18 more variables: SpawnYear <dbl>, kid_year <dbl>, kid_hatchery <chr>,
    ## #   ma_hatchery <chr>, age_spawn <dbl>, mo_dat <dbl>, mo_da <dbl>,
    ## #   program_assn <chr>, Locus <chr>, genecopy_1 <chr>, genecopy_2 <chr>,
    ## #   genotype_number <dbl>, genotype <chr>, typed_final <chr>, sex_f <fct>,
    ## #   year_f <fct>, age_spawn_f <fct>, group <chr>

And now let’s count the different genotypes in the genotype column:

``` r
fg_locus_nest$data[[1]] %>%
  count(genotype)
```

    ## # A tibble: 2 × 2
    ##   genotype     n
    ##   <chr>    <int>
    ## 1 C/C        642
    ## 2 C/T         97

OK, we see that sometimes both homozygotes will not be seen, buy we kind
of need them for making the design matrices. So, what we need is a
function that finds the two alleles, sorts them, and then makes the
genotypes out of them. Let’s just play around with something for doing
that by naming a variable `x` that is `fg_locus_nest$data[[1]]`:

``` r
x <- fg_locus_nest$data[[3]]

# we have to drop missing genotypes:
y <- x$genotype[!is.na(x$genotype)]

# here is some code that gets the unique alleles in sorted order
a <- str_split(y, pattern = "/") %>%
  flatten() %>% 
  as.character() %>%
  unique() %>%
  sort()

a
```

    ## [1] "G"

``` r
# and here is code that makes the sorted possible genotypes vector from uniq_alle
if(length(a) > 1) {
  paste(
    c(a[1], a[1], a[2]), 
    c(a[1], a[2], a[2]), 
    sep = "/"
  )
}
```

So, here is a function that will return the possible genotypes from a
tibble with the genotype column:

``` r
possible_genos <- function(x) {
  y <- x$genotype[!is.na(x$genotype)]
  
  # here is some code that gets the unique alleles in sorted order
  a <- str_split(y, pattern = "/") %>%
    flatten() %>% 
    as.character() %>%
    unique() %>%
    sort()
  
  
  if(length(a) == 2) {
    ret <- paste(
      c(a[1], a[1], a[2]), 
      c(a[1], a[2], a[2]), 
      sep = "/"
    )
  } else if(length(a) == 1) {
    ret <- paste(a[1], a[1], sep = "/")
  } else {
    stop("Locus with 0 or >2 alleles")
  }
  
  ret
}
```

We can test that out like this:

``` r
possible_genos(fg_locus_nest$data[[1]])
```

    ## [1] "C/C" "C/T" "T/T"

And for a mononorphic locus, like this:

``` r
possible_genos(fg_locus_nest$data[[3]])
```

    ## [1] "G/G"

That seems to be working. Note that this function returns a vector of
either length 1 or 3. So, if we want to make a column out of the results
when we apply it to every row of `fg_locus_nest`, then we will have to
make that column a list column to be able to hold it. That is a job for
the `map()` function from the ‘purrr’ package (which is also part of the
tidyverse).

`map()` will apply a function (the `.f` argument) to each row of a
column (the `.x` argument) of a table, and it will return a list column.
So, we can use it inside `mutate()` like this:

``` r
locus_genos <- fg_locus_nest %>%
  mutate(geno_vec = map(.x = data, .f = possible_genos))

locus_genos
```

    ## # A tibble: 24 × 3
    ## # Groups:   Locus_new [24]
    ##    Locus_new                      data                geno_vec 
    ##    <chr>                          <list>              <list>   
    ##  1 Omy_Ch28_greb1_mhap1_position1 <tibble [739 × 27]> <chr [3]>
    ##  2 Omy_Ch28_11658853_F_position1  <tibble [745 × 27]> <chr [1]>
    ##  3 Omy_Ch28_11622532_B_position1  <tibble [746 × 27]> <chr [1]>
    ##  4 Omy_Ch28_11623412_B_position1  <tibble [745 × 27]> <chr [1]>
    ##  5 Omy_Ch28_11641623_D_position1  <tibble [744 × 27]> <chr [1]>
    ##  6 Omy_Ch28_11629294_D_position1  <tibble [732 × 27]> <chr [1]>
    ##  7 Omy_Ch28_greb1_mhap2_position1 <tibble [735 × 27]> <chr [3]>
    ##  8 Omy_Ch28_11625241_C_position1  <tibble [741 × 27]> <chr [3]>
    ##  9 Omy_Ch28_greb1_mhap6_position1 <tibble [745 × 27]> <chr [3]>
    ## 10 Omy_Ch28_greb1_mhap8_position1 <tibble [739 × 27]> <chr [3]>
    ## # ℹ 14 more rows

## *Typed* forms of `map()`

Sometimes you know that your function will return an atomic vector,
i.e. a vector in which each component is a single simple element, and
you probably will know what type they will all be. In that case, there
are a variety of *typed* `map()` functions that will return an atomic
vector (not a list!) and it will check to make sure that they type of
the vector is correct. These different forms are:

- `map_lgl()`: expects to reduce output to a logical vector
- `map_int()`: expects to reduce output to an integer vector
- `map_dbl()`: expects to reduce output to a numeric vector
- `map_chr()`: expects to reduce output to a character vector
- `map_vec()`: expects to reduce output to an atomic vector of any type

### Number of genotypes using `map_int()`

We will demonstrate how to use `map_int()` to make a new column that
gives the number of possible genotypes at each locus.

``` r
locus_genos %>%
  mutate(num_genos = map_int(.x = geno_vec, .f = length))
```

    ## # A tibble: 24 × 4
    ## # Groups:   Locus_new [24]
    ##    Locus_new                      data                geno_vec  num_genos
    ##    <chr>                          <list>              <list>        <int>
    ##  1 Omy_Ch28_greb1_mhap1_position1 <tibble [739 × 27]> <chr [3]>         3
    ##  2 Omy_Ch28_11658853_F_position1  <tibble [745 × 27]> <chr [1]>         1
    ##  3 Omy_Ch28_11622532_B_position1  <tibble [746 × 27]> <chr [1]>         1
    ##  4 Omy_Ch28_11623412_B_position1  <tibble [745 × 27]> <chr [1]>         1
    ##  5 Omy_Ch28_11641623_D_position1  <tibble [744 × 27]> <chr [1]>         1
    ##  6 Omy_Ch28_11629294_D_position1  <tibble [732 × 27]> <chr [1]>         1
    ##  7 Omy_Ch28_greb1_mhap2_position1 <tibble [735 × 27]> <chr [3]>         3
    ##  8 Omy_Ch28_11625241_C_position1  <tibble [741 × 27]> <chr [3]>         3
    ##  9 Omy_Ch28_greb1_mhap6_position1 <tibble [745 × 27]> <chr [3]>         3
    ## 10 Omy_Ch28_greb1_mhap8_position1 <tibble [739 × 27]> <chr [3]>         3
    ## # ℹ 14 more rows

That is cool. Note that some of the loci are monomorphic. We actually
don’t want to hassle with them, because there is nothing interesting to
be found there. So, let’s just filter those out. Note that `map_int()`
can be used anywhere an atomic vector is expected, such as in a
`filter()` statement.

``` r
polymorph_genos <- locus_genos %>%
  filter(map_int(geno_vec, length) == 3) %>%
  select(Locus_new, geno_vec)  # just keep the columns we need later
polymorph_genos
```

    ## # A tibble: 18 × 2
    ## # Groups:   Locus_new [18]
    ##    Locus_new                       geno_vec 
    ##    <chr>                           <list>   
    ##  1 Omy_Ch28_greb1_mhap1_position1  <chr [3]>
    ##  2 Omy_Ch28_greb1_mhap2_position1  <chr [3]>
    ##  3 Omy_Ch28_11625241_C_position1   <chr [3]>
    ##  4 Omy_Ch28_greb1_mhap6_position1  <chr [3]>
    ##  5 Omy_Ch28_greb1_mhap8_position1  <chr [3]>
    ##  6 Omy_Ch28_11678603_G_position1   <chr [3]>
    ##  7 Omy_Ch28_greb1_mhap10_position1 <chr [3]>
    ##  8 Omy_Ch28_11667578_F_position1   <chr [3]>
    ##  9 Omy_Ch28_11683310_H_position1   <chr [3]>
    ## 10 Omy_Ch28_11684744_H_position1   <chr [3]>
    ## 11 Omy_Ch28_11676622_G_position1   <chr [3]>
    ## 12 Omy_Ch28_greb1_mhap2_position2  <chr [3]>
    ## 13 Omy_Ch28_greb1_mhap8_position2  <chr [3]>
    ## 14 Omy_Ch28_greb1_mhap10_position2 <chr [3]>
    ## 15 Omy_Ch28_11683310_H_position2   <chr [3]>
    ## 16 Omy_Ch28_11684744_H_position2   <chr [3]>
    ## 17 Omy_Ch28_11684744_H_position3   <chr [3]>
    ## 18 Omy_Ch28_11684744_H_position4   <chr [3]>

So, that is 18 SNPs.

Ultimately, we will want to run linear models on each of those 18 SNPs
in the two different groups (CV and NH), and when we do that, we need to
have those `geno_vec`’s. So, let’s make a tibble that has all the
tibbles we need for `lmer()` and also the geno_vec. We can toss the
monomorphic loci in the process in one fell swoop by doing an
`inner_join()` that only keeps rows with matching keys:

``` r
dat <- fg_nests %>%
  inner_join(polymorph_genos, by = join_by(Locus_new))
dat
```

    ## # A tibble: 36 × 4
    ## # Groups:   group, Locus_new [36]
    ##    Locus_new                       group data                geno_vec 
    ##    <chr>                           <chr> <list>              <list>   
    ##  1 Omy_Ch28_greb1_mhap1_position1  CV    <tibble [631 × 26]> <chr [3]>
    ##  2 Omy_Ch28_greb1_mhap2_position1  CV    <tibble [627 × 26]> <chr [3]>
    ##  3 Omy_Ch28_11625241_C_position1   CV    <tibble [635 × 26]> <chr [3]>
    ##  4 Omy_Ch28_greb1_mhap6_position1  CV    <tibble [637 × 26]> <chr [3]>
    ##  5 Omy_Ch28_greb1_mhap8_position1  CV    <tibble [635 × 26]> <chr [3]>
    ##  6 Omy_Ch28_11678603_G_position1   CV    <tibble [638 × 26]> <chr [3]>
    ##  7 Omy_Ch28_greb1_mhap10_position1 CV    <tibble [638 × 26]> <chr [3]>
    ##  8 Omy_Ch28_11667578_F_position1   CV    <tibble [637 × 26]> <chr [3]>
    ##  9 Omy_Ch28_11683310_H_position1   CV    <tibble [596 × 26]> <chr [3]>
    ## 10 Omy_Ch28_11684744_H_position1   CV    <tibble [619 × 26]> <chr [3]>
    ## # ℹ 26 more rows

Note that is 36 rows. Just like it should be.

## Dealing with three different models

We know that we want to run three different `lmer()` models on each data
set:

- year_only model: `mo_da ~ I(d_add) + I(d_dom_with_s) + (1|year)`
- sex model: `mo_da ~ I(d_add) + I(d_dom_with_s) + sex + (1|year)`
- age model: `mo_da ~ I(d_add) + I(d_dom_with_s) + age_spawn + (1|year)`

How can we do this? Well there are several ways of doing it, but
probably the tidiest is going to be to lengthen our `dat` tibble to have
three rows for each combination of `group` and `Locus_new`. (This is not
super space-efficient, but it is tidy, and the data sets are not huge).
One way to do this sort of expansion of our data sets is using
`expand_grid()`. To use it, we first make a tibble that has the models
in it. Recall that we made `_f` versions of the variables that are
factors.

``` r
model_tib <- tibble(
  model_name = c("vanilla", "sex", "age"),
  model_formula = c(
    mo_da ~ I(d_add) + I(d_dom_with_s) + (1|year_f),
    mo_da ~ I(d_add) + I(d_dom_with_s) + sex_f + (1|year_f),
    mo_da ~ I(d_add) + I(d_dom_with_s) + age_spawn_f + (1|year_f)
  )
)
model_tib
```

    ## # A tibble: 3 × 2
    ##   model_name model_formula
    ##   <chr>      <list>       
    ## 1 vanilla    <formula>    
    ## 2 sex        <formula>    
    ## 3 age        <formula>

Note that the vector for `model_formula` is a list of model formulas.
`tibble()` is smart enough to know that such a list should be a column
in a tibble. The same cannot be said about `data.frame()`, which is one
reason `tibble()` is superior for dealing with list columns.

Now we can use `expand_grid()` to get us all combinations of the rows in
`dat` with the rows in `model_tib`.

``` r
ready_for_lmer <- expand_grid(dat, model_tib)

ready_for_lmer
```

    ## # A tibble: 108 × 6
    ##    Locus_new                    group data     geno_vec model_name model_formula
    ##    <chr>                        <chr> <list>   <list>   <chr>      <list>       
    ##  1 Omy_Ch28_greb1_mhap1_positi… CV    <tibble> <chr>    vanilla    <formula>    
    ##  2 Omy_Ch28_greb1_mhap1_positi… CV    <tibble> <chr>    sex        <formula>    
    ##  3 Omy_Ch28_greb1_mhap1_positi… CV    <tibble> <chr>    age        <formula>    
    ##  4 Omy_Ch28_greb1_mhap2_positi… CV    <tibble> <chr>    vanilla    <formula>    
    ##  5 Omy_Ch28_greb1_mhap2_positi… CV    <tibble> <chr>    sex        <formula>    
    ##  6 Omy_Ch28_greb1_mhap2_positi… CV    <tibble> <chr>    age        <formula>    
    ##  7 Omy_Ch28_11625241_C_positio… CV    <tibble> <chr>    vanilla    <formula>    
    ##  8 Omy_Ch28_11625241_C_positio… CV    <tibble> <chr>    sex        <formula>    
    ##  9 Omy_Ch28_11625241_C_positio… CV    <tibble> <chr>    age        <formula>    
    ## 10 Omy_Ch28_greb1_mhap6_positi… CV    <tibble> <chr>    vanilla    <formula>    
    ## # ℹ 98 more rows

This is pretty cool. We now have a tibble in which each row contains all
the things we need to run our linear mixed model for certain combination
of Locus_new and group:

- the data
- the vector of possible genotypes
- the model formula

and we have the model name in there to make it easy to see what model it
is.

So, if we had a function that ran `lmer()` with those three inputs and
we could apply that function to each row of `ready_for_lmer`, we could
do that, and put the model results into a new list column.

Since the function for doing that will take three inputs, we have to use
a generalized version of `map()` called `pmap()` (short for “parallel
map”) with which you can pass an arbitrary number of columns—not just 1
(or 2 in the case of `map2()`)—to a function.

Next up is writing a function for running `lmer()`.

## A function to run `lmer()`

This is going to be a function of the data tibble, `D`, the vector of
possible genotypes, `g` and the model formula, `m`. Let’s make some
example variables for testing and stuff, by just grabbing elements from
the first row of our `ready_for_lmer` tibble.

``` r
D <- ready_for_lmer$data[[1]]
g <- ready_for_lmer$geno_vec[[1]]
m <- ready_for_lmer$model_formula[[1]]

D
```

    ## # A tibble: 631 × 26
    ##    indiv   spawner_group sex    year hatchery ArchiveID LENGTH Method collection
    ##    <chr>   <chr>         <chr> <dbl> <chr>    <chr>      <dbl> <chr>  <chr>     
    ##  1 M096525 8/24/16       Male   2017 fyke     M096525      464 FYKE … CH        
    ##  2 M096526 9/1/16        Fema…  2017 fyke     M096526      448 FYKE … FRH       
    ##  3 M096527 9/7/16        Fema…  2017 fyke     M096527      452 FYKE … CH        
    ##  4 M096529 9/11/16       Male   2017 fyke     M096529      445 FYKE … MRH       
    ##  5 M096530 9/13/16       Male   2017 fyke     M096530      380 FYKE … CH        
    ##  6 M096532 9/22/16       Fema…  2017 fyke     M096532      431 FYKE … CH        
    ##  7 M096533 9/25/16       Male   2017 fyke     M096533      480 FYKE … CH        
    ##  8 M096534 9/26/16       Male   2017 fyke     M096534      478 FYKE … CH        
    ##  9 M096535 10/1/16       Male   2017 fyke     M096535      500 FYKE … CH        
    ## 10 M096536 11/30/16      Male   2017 fyke     M096536      500 FYKE … CH        
    ## # ℹ 621 more rows
    ## # ℹ 17 more variables: SpawnYear <dbl>, kid_year <dbl>, kid_hatchery <chr>,
    ## #   ma_hatchery <chr>, age_spawn <dbl>, mo_dat <dbl>, mo_da <dbl>,
    ## #   program_assn <chr>, Locus <chr>, genecopy_1 <chr>, genecopy_2 <chr>,
    ## #   genotype_number <dbl>, genotype <chr>, typed_final <chr>, sex_f <fct>,
    ## #   year_f <fct>, age_spawn_f <fct>

``` r
g
```

    ## [1] "C/C" "C/T" "T/T"

``` r
m
```

    ## mo_da ~ I(d_add) + I(d_dom_with_s) + (1 | year_f)
    ## <environment: 0x15e9c8588>

You don’t have to do this, but I find it really helpful when making a
function to have some test variables.

Looking back at Laura’s notebook we see that the steps in running the
model would be…

1.  Make the design matrix

``` r
dm <- tibble(
  genotype = g,
  ss_ref = c(1, 1, 1),
  d_add = c(0, 1, 2),
  d_dom_with_s = c(0, 1, 0)
)

dm
```

    ## # A tibble: 3 × 4
    ##   genotype ss_ref d_add d_dom_with_s
    ##   <chr>     <dbl> <dbl>        <dbl>
    ## 1 C/C           1     0            0
    ## 2 C/T           1     1            1
    ## 3 T/T           1     2            0

2.  Join the design matrix to the data.

``` r
D2 <- D %>%
  left_join(dm, by = join_by(genotype))
```

3.  Run `lmer()` with the appropriate formula:

``` r
lmer(
  formula = m, 
  data = D2, 
  REML = TRUE
)
```

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

    ## Linear mixed model fit by REML ['lmerModLmerTest']
    ## Formula: m
    ##    Data: D2
    ## REML criterion at convergence: 6326.121
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  year_f   (Intercept) 19.78   
    ##  Residual             36.26   
    ## Number of obs: 631, groups:  year_f, 5
    ## Fixed Effects:
    ## (Intercept)     I(d_add)  
    ##      73.072        7.454  
    ## fit warnings:
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

So, those easy steps can be wrapped up into a function easily, like
this:

``` r
run_lmer <- function(D, g, m) {
  
  dm <- tibble(
    genotype = g,
    ss_ref = c(1, 1, 1),
    d_add = c(0, 1, 2),
    d_dom_with_s = c(0, 1, 0)
  )
  
  D2 <- D %>%
    left_join(dm, by = join_by(genotype))
  
  lmer(
    formula = m, 
    data = D2, 
    REML = TRUE
  )
}
```

This function will return an object that is the return type of `lmer()`.

## Using `pmap()` to do `run_lmer()` over all rows of the tibble

`pmap()` takes an argument `.l` which is a list of the columns you want
to apply to the function. If it is a named list, then the names
correspond to the argument names of the function. This is the easiest
and safest way to run things. If `.l` is not a named list, then the
arguments are done positionally, which can be more error prone.

Anyhoo, in our case, the `.l` list would be like this:

``` r
.l = list(D = data, g = geno_vec, m = model_formula)
```

because we want to assign the `data` column to the `D` argument of
`run_lmer()` and so forth.

So, we run it again inside a `mutate()` like this:

``` r
mod_results <- ready_for_lmer %>%
  mutate(
    lmer = pmap(
      .l = list(D = data, g = geno_vec, m = model_formula),
      .f = run_lmer
    )
  )
```

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

    ## boundary (singular) fit: see help('isSingular')

    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

``` r
mod_results
```

    ## # A tibble: 108 × 7
    ##    Locus_new         group data     geno_vec model_name model_formula lmer      
    ##    <chr>             <chr> <list>   <list>   <chr>      <list>        <list>    
    ##  1 Omy_Ch28_greb1_m… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ##  2 Omy_Ch28_greb1_m… CV    <tibble> <chr>    sex        <formula>     <lmrMdLmT>
    ##  3 Omy_Ch28_greb1_m… CV    <tibble> <chr>    age        <formula>     <lmrMdLmT>
    ##  4 Omy_Ch28_greb1_m… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ##  5 Omy_Ch28_greb1_m… CV    <tibble> <chr>    sex        <formula>     <lmrMdLmT>
    ##  6 Omy_Ch28_greb1_m… CV    <tibble> <chr>    age        <formula>     <lmrMdLmT>
    ##  7 Omy_Ch28_1162524… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ##  8 Omy_Ch28_1162524… CV    <tibble> <chr>    sex        <formula>     <lmrMdLmT>
    ##  9 Omy_Ch28_1162524… CV    <tibble> <chr>    age        <formula>     <lmrMdLmT>
    ## 10 Omy_Ch28_greb1_m… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ## # ℹ 98 more rows

Great! Now we have just run 108 models. How do we go about getting at
the results? We will talk about that later. We can use `map()`-like
functions to get at things, or we can put all the coefficients and
results into tibbles using the ‘broom.mixed’ package.

## One thing you might to see is the summary of all the models

``` r
mrs <- mod_results %>%
  mutate(summary = map(.x = lmer, .f = summary))
mrs
```

    ## # A tibble: 108 × 8
    ##    Locus_new         group data     geno_vec model_name model_formula lmer      
    ##    <chr>             <chr> <list>   <list>   <chr>      <list>        <list>    
    ##  1 Omy_Ch28_greb1_m… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ##  2 Omy_Ch28_greb1_m… CV    <tibble> <chr>    sex        <formula>     <lmrMdLmT>
    ##  3 Omy_Ch28_greb1_m… CV    <tibble> <chr>    age        <formula>     <lmrMdLmT>
    ##  4 Omy_Ch28_greb1_m… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ##  5 Omy_Ch28_greb1_m… CV    <tibble> <chr>    sex        <formula>     <lmrMdLmT>
    ##  6 Omy_Ch28_greb1_m… CV    <tibble> <chr>    age        <formula>     <lmrMdLmT>
    ##  7 Omy_Ch28_1162524… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ##  8 Omy_Ch28_1162524… CV    <tibble> <chr>    sex        <formula>     <lmrMdLmT>
    ##  9 Omy_Ch28_1162524… CV    <tibble> <chr>    age        <formula>     <lmrMdLmT>
    ## 10 Omy_Ch28_greb1_m… CV    <tibble> <chr>    vanilla    <formula>     <lmrMdLmT>
    ## # ℹ 98 more rows
    ## # ℹ 1 more variable: summary <list>

## Or use broom.mixed

``` r
library(broom.mixed)
```

What do the three verbs do?

``` r
x <- mod_results$lmer[[1]]
tidy(x)
```

    ## # A tibble: 4 × 8
    ##   effect   group    term            estimate std.error statistic     df  p.value
    ##   <chr>    <chr>    <chr>              <dbl>     <dbl>     <dbl>  <dbl>    <dbl>
    ## 1 fixed    <NA>     (Intercept)        73.1       9.08      8.05   3.35  0.00265
    ## 2 fixed    <NA>     I(d_add)            7.45      4.80      1.55 625.    0.121  
    ## 3 ran_pars year_f   sd__(Intercept)    19.8      NA        NA     NA    NA      
    ## 4 ran_pars Residual sd__Observation    36.3      NA        NA     NA    NA

``` r
augment(x)
```

    ## # A tibble: 631 × 15
    ##    mo_da `I(d_add)` `I(d_dom_with_s)` year_f .fitted .resid   .hat .cooksd
    ##    <dbl>   <I<dbl>>          <I<dbl>> <fct>    <dbl>  <dbl>  <dbl>   <dbl>
    ##  1    38          0                 0 2017      106.  -68.0 0.0484  0.0940
    ##  2    46          0                 0 2017      106.  -60.0 0.0484  0.0732
    ##  3    52          0                 0 2017      106.  -54.0 0.0484  0.0593
    ##  4    56          0                 0 2017      106.  -50.0 0.0484  0.0508
    ##  5    58          0                 0 2017      106.  -48.0 0.0484  0.0468
    ##  6    67          0                 0 2017      106.  -39.0 0.0484  0.0309
    ##  7    70          0                 0 2017      106.  -36.0 0.0484  0.0264
    ##  8    71          0                 0 2017      106.  -35.0 0.0484  0.0249
    ##  9    76          0                 0 2017      106.  -30.0 0.0484  0.0183
    ## 10   136          0                 0 2017      106.   30.0 0.0484  0.0183
    ## # ℹ 621 more rows
    ## # ℹ 7 more variables: .fixed <dbl>, .mu <dbl>, .offset <dbl>, .sqrtXwt <dbl>,
    ## #   .sqrtrwt <dbl>, .weights <dbl>, .wtres <dbl>

``` r
glance(x)
```

    ## # A tibble: 1 × 7
    ##    nobs sigma logLik   AIC   BIC REMLcrit df.residual
    ##   <int> <dbl>  <dbl> <dbl> <dbl>    <dbl>       <int>
    ## 1   631  36.3 -3163. 6334. 6352.    6326.         627

## Let’s do a quick look at pvalues over loci and groups

``` r
tidied_results <- mod_results %>%
  mutate(tidy = map(.x = lmer, .f = tidy))
```

To access these in a super tidy way, we can unnest them. Let’s do it!

``` r
unnested_tidies <- tidied_results %>%
  select(Locus_new, group, model_name, tidy) %>%
  rename(fish_group = group) %>%
  unnest(tidy)
```

``` r
# then we can look at all of those
unnested_tidies %>%
  filter(term == "I(d_add)") %>%
  arrange(p.value)
```

    ## # A tibble: 108 × 11
    ##    Locus_new         fish_group model_name effect group term  estimate std.error
    ##    <chr>             <chr>      <chr>      <chr>  <chr> <chr>    <dbl>     <dbl>
    ##  1 Omy_Ch28_1168474… NH         sex        fixed  <NA>  I(d_…     89.7      8.78
    ##  2 Omy_Ch28_1168474… NH         sex        fixed  <NA>  I(d_…     89.7      8.78
    ##  3 Omy_Ch28_1168474… NH         sex        fixed  <NA>  I(d_…    -89.7      8.78
    ##  4 Omy_Ch28_1168474… NH         vanilla    fixed  <NA>  I(d_…     88.7      9.40
    ##  5 Omy_Ch28_1168474… NH         vanilla    fixed  <NA>  I(d_…     88.7      9.40
    ##  6 Omy_Ch28_1168474… NH         vanilla    fixed  <NA>  I(d_…    -88.7      9.40
    ##  7 Omy_Ch28_1168331… NH         sex        fixed  <NA>  I(d_…     83.4      9.07
    ##  8 Omy_Ch28_1168331… NH         vanilla    fixed  <NA>  I(d_…     82.9      9.48
    ##  9 Omy_Ch28_1166757… NH         sex        fixed  <NA>  I(d_…    -90.6     10.5 
    ## 10 Omy_Ch28_greb1_m… NH         sex        fixed  <NA>  I(d_…     90.6     10.5 
    ## # ℹ 98 more rows
    ## # ℹ 3 more variables: statistic <dbl>, df <dbl>, p.value <dbl>
