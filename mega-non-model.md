Setting up a run on mega-non-model
================

- [Preliminaries](#preliminaries)
- [Introduction](#introduction)
- [A note on Snakemake](#a-note-on-snakemake)
- [Setting up the units file](#setting-up-the-units-file)
  - [Background on the units file](#background-on-the-units-file)
  - [Making the units file from the dummy fastqs in this tutorial
    repo](#making-the-units-file-from-the-dummy-fastqs-in-this-tutorial-repo)
- [Question for Cassie: What are the adapter
  seqences?](#question-for-cassie-what-are-the-adapter-seqences)
- [Doing Preliminary QC](#doing-preliminary-qc)
  - [Step 1: Make sure I have latest mega-non-model
    code](#step-1-make-sure-i-have-latest-mega-non-model-code)
  - [Step 2. Create a units file and record the parent data
    dir](#step-2-create-a-units-file-and-record-the-parent-data-dir)
  - [Step 3. Do a prelim_qc run with
    mega-non-model](#step-3-do-a-prelim_qc-run-with-mega-non-model)

## Preliminaries

If you want to be able to follow along with some of the R code here it
could be helpful to:

1.  clone this repo:
    `git clone https://github.com/eriqande/mega-bioinf-tutorials.git`.
    This is not essential, but will let you run some functions on the
    dummy data.

If you want to do the snakemake related things you will need to

2.  Have snakemake up and running.
3.  Have a `mega-non-model-wgs-snakeflow` repo somewhere (like on your
    cluster) and have the latest changes in it:

``` sh
# on the main branch of eric's repo (not your own fork) do:
git pull
```

4.  Have a directory with the fastqs and a SampleSheet.csv from a
    NextSeq run.

## Introduction

This is a tutorial intended for my esteemed labmates and members of the
weekly MEGA-bioinformatics group. We are starting to peel a lot of data
off of our own NextSeq machine, and it is important that everyone in the
lab get ready for that firehose.

One of the first orders of business is going to be setting up a simple
way to do some initial/preliminary QA/QC to assess how well these
NextSeq runs are working. Much of what I will be describing related to
this involves things that I recently added to the mega-non-model
workflow to accommodate some of the issues that we foresee running into.
These include:

- Being able to specify a `data_parent_dir` in the config that will be
  prefixed to every path to a fastq file in the units file. This makes
  it easier to store all the necessary files outside of the
  mega-non-model directory. (Note that similar things can be achieved by
  just using symbolic links within a data directory, but we also offer
  the `data_parent_dir` as an option.)
- The addition of a number of new *destination rules*. These are simply
  Snakemake rules that have as input a number of files that one might
  want to produce for a specific purpose without necessarily running the
  whole workflow. We now have new destination rules:
  - `dest_prelim_qc`: for preliminary QA/QC that can be done without any
    read mapping (i.e., running fastp and then summarizing the results
    in various ways).
  - `dest_gvcfs`: for making the gvcfs for every sample, but stopping at
    that point so that other gvcfs can later be made and added to these
    before loading the genomics databases.

We also demonstrate a simple function that we have stored in this
repository (in the `R` directory) that is tailored for converting our
NextSeq sample sheets and sequencing outputs into a `units.tsv` file for
running mega-non-model.

## A note on Snakemake

The mamba package manager has undergone a few breaking changes, so the
newest versions of Snakemake have been using the latest versions of
conda (which use libmamba under the hood now.) I set that up on one of
my laptops, but I am currently still just using Snakemake 8.5.3 on SEDNA
and on my work laptop, and that seems to be working pretty well with the
old version of mamba that I have.

However, we might want to be prepared for glitches for new users.

## Setting up the units file

We are going to focus on these two first, because those are the ones
that are needed to do preliminary qc (a chromosomes.tsv and
scaffold_groups.tsv file will need to be there, but won’t be relevant.)

### Background on the units file

Making the units file is a critical step, and it is easily done in R
from the file listing of your fastqs, and from the SampleSheet.csv that
comes off the machine. Here we just do a little tour of the background.

First, let’s review the mega-non-model instructions regarding this file.
See it
[here](https://github.com/eriqande/mega-non-model-wgs-snakeflow?tab=readme-ov-file#what-the-user-must-do-and-values-to-be-set-etc).

To make this file, I got most of the necessary ingredients by using `ls`
on SEDNA, where I have stored the files:

``` sh
[sedna: eriq]--% pwd
/share/swfsc/eriq
[sedna: eriq]--% ls -l Landscape_Genomics/*/*.fastq.gz
-rw-rw-r-- 1 eanderson swfsc 1673418258 Mar  2 14:11 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004147_M028_1A_S62_R1_001.fastq.gz
-rw-rw-r-- 1 eanderson swfsc 1779605042 Mar  2 14:13 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004147_M028_1A_S62_R2_001.fastq.gz
-rw-rw-r-- 1 eanderson swfsc 1608758499 Mar  2 14:11 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004182_M028_5D_S61_R1_001.fastq.gz
-rw-rw-r-- 1 eanderson swfsc 1698883931 Mar  2 14:13 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004182_M028_5D_S61_R2_001.fastq.gz
-rw-rw-r-- 1 eanderson swfsc  845174892 Mar  2 14:11 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004225_M028_10G_S74_R1_001.fastq.gz
-rw-rw-r-- 1 eanderson swfsc  884389611 Mar  2 14:12 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004225_M028_10G_S74_R2_001.fastq.gz
-rw-rw-r-- 1 eanderson swfsc  799256079 Mar  2 14:11 Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004226_M028_10H_S75_R1_001.fastq.gz
...
```

and so forth.

The SampleSheet.csv that should be in the same directory as the fastqs
looks [like
this](https://raw.githubusercontent.com/eriqande/mega-bioinf-tutorials/refs/heads/main/data/mega-non/SampleSheet.csv)

Ultimately, we want to combine the information in the file listing with
what is in the sample sheet to get columns that have names like:

    sample  unit    library flowcell    platform    lane    sample_id   barcode fq1 fq2 kb1 kb2

Doing this involves just a little rigamorale in R. I have wrapped up all
the steps in some functions that are stored at:
<https://raw.githubusercontent.com/eriqande/mega-bioinf-tutorials/refs/heads/main/R/mega-non-setup.R>.

Let’s have a quick look at those.

### Making the units file from the dummy fastqs in this tutorial repo

I made a bunch of fastq files with no content and saved them into this
repo in the directory:
`"data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/`.
I also put a `SampleSheet.csv` in that directory.

I tossed out a few of the fastq file pairs, and I also deleted a few of
the samples from the `SampleSheet.csv` in the directory so that everyone
can see what happens when there are fastqs that are not in the
SampleSheet, or samples in the SampleSheet that don’t have fastqs.
(There are some warnings issued and it returns a report).

The main function for this is `write_units_file()` and it takes the
relative path to the directory holding all the fastqs. (As we will see
later, it makes good sense to put all of these different run directories
within a single parent data directory). Let’s run it!

``` r
functions <- "https://raw.githubusercontent.com/eriqande/mega-bioinf-tutorials/refs/heads/main/R/mega-non-setup.R"
source(functions)  # this defines the needed functions
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# then call it
dump <- write_units_file("data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX")
```

    ## Warning in
    ## write_units_file("data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX"):
    ## Not all samples in SampleSheet.csv have fastq paths.  Check return tibble
    ## $missing_paths for info.

    ## Warning in
    ## write_units_file("data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX"):
    ## Not all fastq paths have corresponding entries in SampleSheet.csv.  Check
    ## return tibble $missing_barcodes for info.

We see that it gave some warnings. But it also created a units file that
is inside the fastq directory. This is in:
`data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/units.tsv`,
and it looks like:

``` sh
head data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/units.tsv
```

    ## sample   unit    library flowcell    platform    lane    sample_id   barcode fq1 fq2 kb1 kb2
    ## M004147_M028_1A  1   M004147_M028_1A_TAGCTGGC_TACCGGCT   LGRun1  ILLUMINA    1   M004147_M028_1A TAGCTGGC+TACCGGCT   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004147_M028_1A_S62_R1_001.fastq.gz data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004147_M028_1A_S62_R2_001.fastq.gz 0   0
    ## M004182_M028_5D  1   M004182_M028_5D_GTTGGCGT_CAGCAGTC   LGRun1  ILLUMINA    1   M004182_M028_5D GTTGGCGT+CAGCAGTC   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004182_M028_5D_S61_R1_001.fastq.gz data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004182_M028_5D_S61_R2_001.fastq.gz 0   0
    ## M004225_M028_10G 1   M004225_M028_10G_ATTGGACA_TTGCTGGA  LGRun1  ILLUMINA    1   M004225_M028_10G    ATTGGACA+TTGCTGGA   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004225_M028_10G_S74_R1_001.fastq.gz    data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004225_M028_10G_S74_R2_001.fastq.gz    0   0
    ## M004226_M028_10H 1   M004226_M028_10H_ACCCGTTG_CAGCTTCG  LGRun1  ILLUMINA    1   M004226_M028_10H    ACCCGTTG+CAGCTTCG   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004226_M028_10H_S75_R1_001.fastq.gz    data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004226_M028_10H_S75_R2_001.fastq.gz    0   0
    ## M004444_M031_2B  1   M004444_M031_2B_CCTTCCAT_CACGCAAT   LGRun1  ILLUMINA    1   M004444_M031_2B CCTTCCAT+CACGCAAT   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004444_M031_2B_S30_R1_001.fastq.gz data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004444_M031_2B_S30_R2_001.fastq.gz 0   0
    ## M004446_M031_2D  1   M004446_M031_2D_TTCCAGGT_CAGTGCTT   LGRun1  ILLUMINA    1   M004446_M031_2D TTCCAGGT+CAGTGCTT   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004446_M031_2D_S3_R1_001.fastq.gz  data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004446_M031_2D_S3_R2_001.fastq.gz  0   0
    ## M004452_M031_3B  1   M004452_M031_3B_CCAGTATC_AGCGAGAT   LGRun1  ILLUMINA    1   M004452_M031_3B CCAGTATC+AGCGAGAT   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004452_M031_3B_S10_R1_001.fastq.gz data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004452_M031_3B_S10_R2_001.fastq.gz 0   0
    ## M004453_M031_3C  1   M004453_M031_3C_GTCAGTCA_CGGCATTA   LGRun1  ILLUMINA    1   M004453_M031_3C GTCAGTCA+CGGCATTA   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004453_M031_3C_S29_R1_001.fastq.gz data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004453_M031_3C_S29_R2_001.fastq.gz 0   0
    ## M004465_M031_4G  1   M004465_M031_4G_ATCTGACC_CTCGACTT   LGRun1  ILLUMINA    1   M004465_M031_4G ATCTGACC+CTCGACTT   data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004465_M031_4G_S23_R1_001.fastq.gz data/mega-non/empty_files/Landscape_Genomics/250301_VH02170_2_2227GLJNX/M004465_M031_4G_S23_R2_001.fastq.gz 0   0

If you want to investigate the errors, you can see the SampleSheet
entries that don’t have fastqs with this:

``` r
dump$missing_paths
```

    ## # A tibble: 2 × 15
    ##   sample         kb1   kb2 fq1   fq2    lane sample_id platform flowcell library
    ##   <chr>        <dbl> <dbl> <chr> <chr> <int> <chr>     <chr>    <chr>    <chr>  
    ## 1 M076895_M77…    NA    NA <NA>  <NA>     NA <NA>      <NA>     LGRun1   M07689…
    ## 2 M076891_M77…    NA    NA <NA>  <NA>     NA <NA>      <NA>     LGRun1   M07689…
    ## # ℹ 5 more variables: LibraryPrepKitName <lgl>, IndexAdapterKitName <lgl>,
    ## #   Index <chr>, Index2 <chr>, barcode <chr>

And you can see the fastqs that don’t have corresponding SampleSheet
entries with this:

``` r
dump$missing_barcodes
```

    ## # A tibble: 4 × 15
    ##   sample         kb1   kb2 fq1   fq2    lane sample_id platform flowcell library
    ##   <chr>        <dbl> <dbl> <chr> <chr> <int> <chr>     <chr>    <chr>    <chr>  
    ## 1 M076712_M77…     0     0 data… data…     1 M076712_… ILLUMINA <NA>     <NA>   
    ## 2 M077119_M77…     0     0 data… data…     1 M077119_… ILLUMINA <NA>     <NA>   
    ## 3 M077120_M77…     0     0 data… data…     1 M077120_… ILLUMINA <NA>     <NA>   
    ## 4 M079477_M80…     0     0 data… data…     1 M079477_… ILLUMINA <NA>     <NA>   
    ## # ℹ 5 more variables: LibraryPrepKitName <lgl>, IndexAdapterKitName <lgl>,
    ## #   Index <chr>, Index2 <chr>, barcode <chr>

## Question for Cassie: What are the adapter seqences?

I just want to check and make sure that we have the right ones in the
workflow. (Though with fastp it doesn’t make too much difference).

## Doing Preliminary QC

I have modified the mega-non-model workflow to have a “destination rule”
that simply requests that fastp get run on everything and then all those
results get displayed with multiqc.

These steps do not require a reference genome or a list of chromosomes
and scaffold groups, etc., so that it is not necessary to customize a
config file for them. You just have to be able to:

1.  tell snakemake where the data folder is
2.  tell snakemake where the units.tsv folder is

Here I am going to demonstrate how this all works in a series of steps
on SEDNA. It should work just fine on hummingbird or even megabox, as
well.

### Step 1: Make sure I have latest mega-non-model code

Just getting my new additions.

``` sh
# note I am in my "test" repo for mega-non-model
(snakemake-8.5.3) [sedna: mega-non-model-wgs-snakeflow]--% pwd
/home/eanderson/scratch/TEST/mega-non-model-wgs-snakeflow

# pull changes into main
(snakemake-8.5.3) [sedna: mega-non-model-wgs-snakeflow]--% git pull
remote: Enumerating objects: 6, done.
remote: Counting objects: 100% (6/6), done.
remote: Total 6 (delta 5), reused 6 (delta 5), pack-reused 0 (from 0)
Unpacking objects: 100% (6/6), done.
From github.com:eriqande/mega-non-model-wgs-snakeflow
   f748c07..a6fe699  main                         -> origin/main
   6b46b7f..a6fe699  add-features-for-lab-nextseq -> origin/add-features-for-lab-nextseq
Updating f748c07..a6fe699
Fast-forward
 .test/config/config.yaml             |  9 +++++++++
 smk-8-slurm/jdb/sedna/config.yaml    |  4 +---
 workflow/rules/common.smk            | 11 ++++++++++-
 workflow/rules/destination-rules.smk |  7 +++++++
 workflow/rules/qc.smk                | 23 +++++++++++++++++++++++
 5 files changed, 50 insertions(+), 4 deletions(-)
```

### Step 2. Create a units file and record the parent data dir

``` sh
# note that I am in what I consider to be my parent data dir. So the output
# here is the absolute path to my parent_data_dir
[node10: Landscape_Genomics]--% pwd
/share/swfsc/eriq/Landscape_Genomics

# do a short tail listing of the run directory we want to make a units file for
[node10: Landscape_Genomics]--% ls 250301_VH02170_2_2227GLJNX/* | tail
250301_VH02170_2_2227GLJNX/M086769_M881_5F_S8_R2_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086773_M881_6B_S7_R1_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086773_M881_6B_S7_R2_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086813_M881_11B_S17_R1_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086813_M881_11B_S17_R2_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086821_M881_12B_S37_R1_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086821_M881_12B_S37_R2_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086826_M881_12G_S44_R1_001.fastq.gz
250301_VH02170_2_2227GLJNX/M086826_M881_12G_S44_R2_001.fastq.gz
250301_VH02170_2_2227GLJNX/SampleSheet.csv
```

Note that the SampleSheet.csv is in there.

Now, we use R to create the units file:

``` r
[node10: Landscape_Genomics]--% module load R/4.0.3
[node10: Landscape_Genomics]--% R

# now we are in R:
functions <- "https://raw.githubusercontent.com/eriqande/mega-bioinf-tutorials/refs/heads/main/R/mega-non-setup.R"
source(functions)  # this defines the needed functions

# then call the function
dump <- write_units_file(fastq_dir = "250301_VH02170_2_2227GLJNX")
```

### Step 3. Do a prelim_qc run with mega-non-model

For this, we are back in our mega-non-model directory:

``` sh
[node10: mega-non-model-wgs-snakeflow]--% pwd
/home/eanderson/scratch/TEST/mega-non-model-wgs-snakeflow
```

We will do a dry-run. We use the default config file in the .test
directory, but we update the units file location and the parent_data_dir
on the command line.

``` sh
[node10: mega-non-model-wgs-snakeflow]--% conda activate snakemake-8.5.3


snakemake -np --cores 20 --use-conda  dest_prelim_qc --configfile .test/config/config.yaml   --config units= /share/swfsc/eriq/Landscape_Genomics/250301_VH02170_2_2227GLJNX/units.tsv  parent_data_dir=/share/swfsc/eriq/Landscape_Genomics/
```

Notice that you need to have a trailing slash on the parent_data_dir
when you specify it.
