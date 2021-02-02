---
tags: ggg, ggg2021, ggg201b
---
# Lab 4 outline - even more variant calling; snakemake wildcards. - GGG 201(b) lab, Jan 29, 2021

## last week outline + additions

last Friday was mostly "what are we doing with variant calling, and this workflow specifically?"

We'll finish that discussion today.

## but first... homework.

I plan to post homework today, to be due Tuesday, February 9th. Does this conflict with other classes?

### last week (week 3)

assumptions and scope:
- which do we want? SNPs, short indels, long indels/structural variation
- accurate short reads work well for SNPs
- inaccurate long reads vs accurate long reads

workflow steps and how they affect results -
 * affected by biology (genome size, structure, etc.)
 * affected by sequencing details (long reads vs short, etc.)
 * affected by parameters (depth, diploidy, etc.)
 * not affected at all


Outline of likely topics:
* how shotgun sequencing works
    * sample prep
    * Illumina + PCR
    * Single molecule sequencing & long reads
* key assumptions underlying our sequence analysis: "quantitative" or "digital" sequencing
    * variant calling
    * assembly
    * RNAseq
    * vs metagenomics
* our goal(s) in variant calling

### week 4

* handling the data
    * quality analysis and error removal
    * alignment
    * organizing the data
    * calling variants via pileup
        * (what other methods could be used?)
    * summarizing variants, and interpreting VCF
    * visualization and "gut check"
* how does biology mess with us?
    * ploidy
    * repetitive sequence
    * heterozygosity/strain variation

## week 4 - more snakemake: wildcards!

Let's make this workflow a bit more generic/less sample specific!

[![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/ngs-docs/2021-ggg-201b-variant-calling/week3?urlpath=rstudio)

If you look at [our week2 workflow](https://github.com/ngs-docs/2021-ggg-201b-variant-calling/blob/week3/Snakefile), you'll see that we have a specific sample names in there (`SRR2584857_1`). How do we make a workflow that runs on more than one sample?

The answer is "wildcards". Snakemake can *automatically* figure out the sample name, given a few hints.

To make use of wildcards, you replace the sample name in the input and output for a rule with "{substitute}".

Let's do this together for map_reads, and then you can do it on your own for sam_to_bam and sort_bam.

```python
rule map_reads:
    conda: "env-minimap.yml"
    input:
        ref="ecoli-rel606.fa.gz",
        reads="{sample}.fastq.gz"
    output: "{sample}.ecoli-rel606.sam"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """
```

A few simple rules for wildcards:

* wildcard variable names ('substitute', above) are specific to a rule.
    * The same variable name must be used within each rule's inputs, outputs, etc.
    * Variable names don't mean anything, so you can use 'sample' or whatever.
    * snakemake tells you what it's substituting when it runs - look at the "wildcards" output.
* you must also always provide rules or filenames when you are using wildcards - snakemake doesn't (can't) automatically figure out what wildcards are appropriate.

Question: how does snakemake figure out what the wildcard values should be?

(The answer is "pattern matching" - it guesses based on what it's being asked to do, but only when it can do so with 100% reliability, under its rules! ...a general feature of computers :)

"Concrete" vs "wildcard" rules -
* concrete rules ask for / produce a specific filename, without wildcards.
* wildcard rules provide a _recipe_ for producing a general filename, using wildcards.

Let's try doing this for call_variants. ...what goes wrong?
