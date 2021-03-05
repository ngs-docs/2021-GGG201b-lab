---
tags: ggg, ggg2021, ggg298
---
# GGG 201(b), Lab Homework 2 - 2021

Due by 7pm, Monday March 15th

Please note that you may freely work together with others, but you must hand it in separately.

Remember to execute all of the below in an srun session!

## Part 1: Add more samples for variant calling to HW #1

Starting in your first HW repository, add the following samples to your variant calling workflow:

`SRR2584857_1.fastq.gz` - `https://osf.io/aksmc/download`
`SRR2584405_1.fastq.gz` - `https://osf.io/7qek6/download`
`SRR2584404_1.fastq.gz` - `https://osf.io/s24ky/download`
`SRR2584403_1.fastq.gz` - `https://osf.io/6jx7z/download`

Your snakemake workflow should download these files if they don't exist, and produce separate .vcf files for each one.

There are multiple ways to do this --

You can do this using [trick #1 in my blog post](http://ivory.idyll.org/blog/2020-snakemake-hacks-collections-files.html) if you like. 

Simpler and easier-to-debug alternatives:
* create four rules that each provide an `output:` of a different one of the raw data files, and provide a wget or curl command in the associated`shell:` block that produces each file.
* create a single rule that provides a single `output:` block that contains all four of the raw data file names, and write four separate commands (one on each line) in the `shell:` block that download the four files.

Side note: you can now run with `-j 4` to run it ~four times faster.

### Outputs

You should end up with four (different) VCF files.

### Things to check

* Make sure that `snakemake --delete-all-output` removes all of the generated files.
* Make sure that you can rerun the pipeline from scratch without any errors.

### Handing it in

Run `git commit -am` and `git push` to submit. Make sure that your latest changes are reflected on your github site.

## Part 2: Where does assembly break down?

After accepting the HW #2 assignment URL (see Canvas or Slack), please clone the new github repository for assembly into your farm account.

Adjust the Snakefile to generate assemblies and quast evaluations for three different subsets (I suggest you keep them between 0 and 2000000 in size).

After sucessfully generating three different assemblies and quast evaluations, fill out [this form](https://docs.google.com/forms/d/e/1FAIpQLScCmsRdHIujMcwYp-vdToe8yd6d58rDLNVmHkiAZ_YR5FjgXA/viewform) with your resulting numbers.

Also please commit and push your Snakefile.
