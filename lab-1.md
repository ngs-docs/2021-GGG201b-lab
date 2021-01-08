---
tags: ggg, ggg2021, ggg201b
---

[toc]

# GGG 201b, jan 2021 - lab 1 - NOTES

## Friday Lab Outline - 1/10

* introductions and a short presentation
* [syllabus](https://hackmd.io/wnAlw5Y6QRu4kfWiri9Cwg?view)
* our first workflow: variant calling!!!!!!

## Day 1: Variant Calling

learning goals:
- introduce basic concepts of variant calling
- understand the basics of running shell commands via snakemake
- download and examine some data

### A very brief discussion of variant calling

(Back to the presentation)

### Start a binder

https://github.com/ngs-docs/2021-ggg-201b-variant-calling

[![Binder](https://mybinder.org/badge_logo.svg)](https://binder.pangeo.io/v2/gh/ngs-docs/2021-ggg-201b-variant-calling/week1?urlpath=rstudio)

wait 20 or 30 seconds...

### Open a terminal

...and set a prompt:
```
PS1='$ '
```

### Run some stuff!

There are lots of "rules" in the Snakefile [(view it on GitHub, too)](https://github.com/ngs-docs/2021-ggg-201b-variant-calling/week1/Snakefile)! Try:

```
grep rule Snakefile
```

and

```
nano -ET4 Snakefile
```
to examine them. To exit type `Control + x`

* what is the structure of a rule?
  * 'rule' block
  * indentation
  * "shell" and "conda" as well as others
* what does shell mean here?
* what is the "conda" block? we'll talk about that in a sec...

----

Let's run one! Try the following command:
```
snakemake -p -j 1 map_reads
```

This says, "run the shell command in the Snakefile under the 'map_reads' rule".

`-p` means "show the command that you're running.

`-j 1`

It will fail! Why?

Because we don't have the minimap software installed.

If we say `--use-conda`, snakemake will automatically install the right
software (specified in `env-minimap.yml`).

Let's try that:

```
snakemake -p -j 1 map_reads --use-conda
```

---

This also fails! Why?

Because we don't have any data! We need to download some stuff and prepare it.

Start with

```
snakemake -p -j 1 --use-conda download_data
```

ahh... cool green!

What does this command do? tl;dr creates the file `SRR2584857_1.fastq.gz`.
(What's in this file?)

Now run some more -- one at a time.

```
snakemake -p -j 1 --use-conda download_genome
```
This creates the file `ecoli-rel606.fa.gz.

Now run:
```
snakemake -p -j 1 --use-conda map_reads
snakemake -p -j 1 --use-conda sam_to_bam
snakemake -p -j 1 --use-conda sort_bam
snakemake -p -j 1 --use-conda call_variants
```

This runs a complete (if fairly simple :) variant calling workflow.

The result is a file called `SRR2584857_1.ecoli-rel606.vcf`.

View the results in the RStudio window by clicking on that file.

You can also go to the shell prompt and execute:

```
samtools index SRR2584857_1.ecoli-rel606.bam.sorted
samtools tview -p ecoli:4314717 --reference ecoli-rel606.fa SRR2584857_1.ecoli-rel606.bam.sorted
```
which will show you the actual aligned reads.

A few things to discuss:

* these are just shell commands, with a bit of "decoration". You could run them yourself if you wanted!
* order of the rules in the Snakefile doesn't matter
* rules can have one or more shell commands
* `snakemake -p` prints the command
* `-j 1` says "use only one CPU"
* `--use-conda` says to install software as specified.
* red if it fails (try breaking one command :)
* it's all case sensitive...
* tabs and spacing matter.
* why are the files named the way they are?
* if you tell snakemake what software a step needs, it will install that software for you (and manage it for you)

Some foreshadowing:

* wouldn't it be nice to not have to run each rule one by one?
* wouldn't it be nice if the command didn't rerun if it had already run

### Upgrading our Snakefile by adding 'output:'

Edit `Snakefile` and add:
and add `output: "SRR2584857_1.fastq.gz"` to the `download_data` rule.

Now try running the rule: `snakemake -p SRR2584857_1.fastq.gz`

Run it again. Hey look, it doesn't do anything! Why??

Try removing the file: `rm SRR2584857_1.fastq.gz`. Now run the rule again.

### Upgrading our Snakefile by adding `input:`.

To the `map_reads` rule, add:

`input: "SRR2584857_1.fastq.gz", "ecoli-rel606.fa"`

What does this do? (And does it work?)

### Rewrite the rule shell blocks to use `{input}` and `{output}`

For each rule where you have defined `input:` and `output:` you can replace the
filenames in the `shell:` block with `{input}` and `{output}` as appropriate.

### End of class recap

* Snakefile contains a snakemake workflow definition
* the rules specify steps in the workflow
* At the moment (and in general), they run shell commands
* You can "decorate" the rules to tell snakemake how they depend on each other.
