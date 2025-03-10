Setting up a run on mega-non-model
================

- [Introduction](#introduction)
- [A note on Snakemake](#a-note-on-snakemake)

## Introduction

This is a tutorial intended for my esteemed labmates and members of the
weekly MEGA-bioinformatics group. We are starting to peel a lot of data
off of our own NextSeq machine, and it is important that everyone in the
lab get ready for that firehose.

One of the first orders of business is going to be setting up a simple
way to do some initial/preliminary QA/QC to assess how well these
NextSeq runs are working. Much of what I will be describing related to
this involve things that I recently added to the mega-non-model workflow
to accommodate some of the issues that we foresee running into. These
include:

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
