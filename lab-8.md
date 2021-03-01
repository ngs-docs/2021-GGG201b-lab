---
tags: ggg, ggg2021, ggg201b
---
# Lab 8 outline - finishing off assembly and variant calling! GGG 201(b) lab, Feb 26, 2021

[toc]

## Note up front

I plan to assign at least one more homework, this afternoon. Due date March 9th.

## Preparation

Log into farm & make a new copy of HW 1 - 

```
git clone https://github.com/ngs-docs/2021-ggg-201b-lab-hw1 ~/lab8-vc
```

Grab a compute node:
```
srun --nodes=1 -p high2 -t 2:00:00 -c 8 --mem 10GB --pty /bin/bash
```

Activate your `vc` conda environment (which you created as part of HW #1; you can create it now with `conda env create -n vc -f ~/lab8-vc/environment.yml`) --

```
conda activate vc
```

You'll need to install a specific version of samtools:
```
conda install -y samtools=1.10
```


Now run the snakemake workflow --

```
cd ~/lab8-vc/
snakemake -j 1 --use-conda
```

Note, this is the [same old workflow](https://github.com/ngs-docs/2021-ggg-201b-lab-hw1/blob/main/Snakefile) from hw1.

## Variant calling and VCF files revisited

### Getting coverage stats

We haven't really talked about the uncovered regions - how do we get at those?

`samtools depth` will report per-base coverage.

```
samtools depth -aa SRR2584857_1.ecoli-rel606.bam.sorted > SRR2584857_1.ecoli-rel606.depth
head SRR2584857_1.ecoli-rel606.depth
```
and these are now accessible to spreadsheets, Python, and R.

### Finding mapped and unmapped reads

samtools supports a bunch of different "flags" for getting at various types of reads - see this [helpful site](https://broadinstitute.github.io/picard/explain-flags.html).

Both the `view` and `bam2fq` commands support probing these flags with `-f` (flag is present) and `-F` (flag is not present).

For example, if you [look at the Snakefile](https://github.com/ngs-docs/2021-ggg-201b-lab-hw1/blob/main/Snakefile#L37), we use `samtools view -F 4` to find all _mapped_ reads (this says we want "SAM entries that do not have the unampped flag" :shrug:)

You can export all of the mapped reads with `bam2fq` like so:
```
samtools bam2fq -F 4 SRR2584857_1.ecoli-rel606.bam
```

How many unmapped reads are there?

(What's going on here?)

### Specific vs sensitive variant calling

If you look at the VCF file,

```
grep -v ^# *.vcf
```

you'll see lines that look like this
```
ecoli   920514  .       T       C       69      .       DP=3;VDB=0.0900131;SGB=-0.511536;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,1,2;MQ=60      GT:PL   1/1:99,9,0
ecoli   4141441 .       C       T       41.4148 .       DP=3;VDB=0.32;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=60  GT:PL   1/1:71,6,0
```
where the fields are chromosome, position, allele name, reference allele, variant allele, quality score, whether they have been filtered, and an INFO field.

What are these "INFO" fields, e.g.
`DP=3;VDB=0.0900131;SGB=-0.511536;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,1,2;MQ=60`?

The VCF fields are specified in the [VCF specs](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf)

Here `DP=3` means depth 3, for example.

You can filter variant calls based on these various elements using bcftools, like so:

```
bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' 
```

* QUAL is mapping quality, aka ["phred score"](https://en.wikipedia.org/wiki/Phred_quality_score) - here we are excluding base calls that are estimated to be less than 99.99% accurate (QUAL < 40)
* DP is depth - here we are excluding depth < 10
* GT is genotype - here we are excluding heterozygous alleles

Let's modify the Snakefile to produce both _sensitive_ and _specific_ variant calls - the sensitive calls with no filtering, the specific calls with the above filtering.

* split the variant calling rule in two
* write a new rule

See [Torsten Seeman's blog post](http://thegenomefactory.blogspot.com/2018/10/a-unix-one-liner-to-call-bacterial.html) for the source of  these parameters.)

### Looking at variants

I wrote a convenient little script to look at variants in the VCF file - [view the source code](https://github.com/ngs-docs/2021-GGG201b-lab/blob/latest/lab-8.run-tview.py), and then download it,
```
curl -L https://raw.githubusercontent.com/ngs-docs/2021-GGG201b-lab/latest/lab-8.run-tview.py > run-tview.py
chmod +x run-tview.py
```

and run it like so:
```
./run-tview.py SRR2584857_1.ecoli-rel606
```

## What didn't assemble?

We now have the tools to look at what reads _didn't_ assemble --
```
mkdir ~/lab8-assembly/
cd ~/lab8-assembly/
ln -fs ~ctbrown/data/ggg201b/SRR2584857_*.fastq.gz ./
cp ~ctbrown/data/ggg201b/SRR2584857.full.contigs.fa .
```

What we need to do is:
* try to map all the reads to the assembly
* then, output the reads that _didn't_ map successfully

:notes: shall we build a Snakefile? :notes:


### snakefile

```
rule map_reads:
    input:
        reads_1 = "SRR2584857_1.fastq.gz",
        reads_2 = "SRR2584857_2.fastq.gz",
        ref = "SRR2584857.full.contigs.fa",
    output:
        "SRR2584857.full.sam",
    conda: "env-minimap.yml",
    shell: """
        minimap2 -ax sr {input.ref} {input.reads_1} {input.reads_2} > {output}
    """
```

### commands:

```
% samtools view -f 4 SRR2584857.full.sam | wc -l
33387
% samtools view -F 4 SRR2584857.full.sam | wc -l
4226686
% python
Python 3.9.2 | packaged by conda-forge | (default, Feb 21 2021, 05:02:46) 
[GCC 9.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> unaligned = 33387
>>> aligned = 4226686
>>> unaligned / (aligned + unaligned)
0.007837189644402807
>>> 
```

so, about 0.78% of the reads did NOT align back to the assembly. Personally, I'm OK with ignoring those 33,387 reads.

## Extra: Exploring interval math

[Bedtools](https://bedtools.readthedocs.io/en/latest/index.html) is a ridiculously useful tool that you can use to do all sorts of genome arithmetic, like "find overlapping features" and "find nearest features".

Let's take a look at an example -- here's the GFF file we produced in lab 7 using prokka:
```
cp ~ctbrown/data/ggg201b/SRR2584857.full.gff .
```

Select out the hsrA entry into its own gff file:
```
grep hsrA SRR2584857.full.gff > hsrA.gff
```

and now get the relevant FASTA:

```
bedtools getfasta -fi SRR2584857.full.contigs.fa -bed hsrA.gff
```

## Extra: exploring samtools view for coordinates

How might we get at the reads that align to this specific feature?

1. use `samtools view` to make a new BAM for the hsrA coordinates.
2. use `samtools bam2fq` on that resulting BAM file.

### commands

```
# convert SAM file into BAM file
samtools view -b SRR2584857.full.sam > SRR2584857.full.bam

# sort BAM file by position of read in assembly
samtools sort SRR2584857.full.bam > SRR2584857.full.bam.sorted

# index the BAM for random retrieval
samtools index SRR2584857.full.bam.sorted

# create a BAM from reads that overlap just the position of hsrA in the assembly
samtools view SRR2584857.full.bam.sorted k119_11:86-1133 -b > hsrA.bam

# and output associated FASTQ records
samtools bam2fq hsrA.bam
```

## Extra: Integrative Genomics Viewer

[IGV](http://software.broadinstitute.org/software/igv/) is a really useful way to view genomes and annotations.

## Extra: An example of some code to probe VCF files

From some COVID work I'm doing:

* depth profile and stats reporting - what's the right way to *summarize* depth!?
* digging into specific variants with Python code
