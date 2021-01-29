---
tags: ggg, ggg2021, ggg201b
---
# Lab 3 outline - variant calling, discussed. - GGG 201(b) lab, Jan 22, 2021

## last week reminder

what we did last week -
- defined an end-to-end workflow that led us from a sample to a set of called variants
- we did this by connecting rules to each other (show Snakefile)

## today: lab 3/discussion and whiteboarding

this Friday will be mostly "what are we doing with variant calling, and this workflow specifically?"

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

See also [hand-drawn notes (PDF)](lab-3-notes.pdf)...
