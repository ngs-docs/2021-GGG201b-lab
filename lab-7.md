---
tags: ggg, ggg2021, ggg201b
---
# Lab 7 outline - assembly, part 3! GGG 201(b) lab, Feb 19, 2021

[toc]

## Preparation

Log into farm.

### install a new [sourmash](https://sourmash.readthedocs.io/en/latest/) environment for k-mer exploration

create a new sourmash env:
```
conda create -n smash -c bioconda -c conda-forge -y sourmash 
```

now, activate your sourmash environment and install upgraded sourmash and termplotlib:
```
conda activate smash
pip install --pre sourmash==4.0.0rc1
pip install termplotlib==0.3.4
```

### Grab a compute note:

```
srun --nodes=1 -p high2 -t 2:00:00 -c 8 --mem 10GB --pty /bin/bash
```

Activate your `assembly` environment, which contains snakemake and a few other things:
```
conda activate assembly
```

### get today's assembly workflow

next, grab [the workflow](https://github.com/ngs-docs/2021-ggg-201b-assembly) from github and put it in ~/lab7 --
```
git clone https://github.com/ngs-docs/2021-ggg-201b-assembly ~/lab7
cd ~/lab7
cp ~ctbrown/data/ggg201b/SRR2584857_*.fastq.gz ./
```

and now start running it:
```
snakemake -j 8 --use-conda
```
this will take a while, because it needs to install a lot of software!

(preparation **stops here**)

## Notes on Snakefile for quast and prokka

(We'll start class **here**.)

the goal of this snakefile is to run [quast (for evaluating assemblies)](http://bioinf.spbau.ru/quast) and [prokka (for annotating them)](https://github.com/tseemann/prokka).

note, to run things, make sure you're in the `assembly` environment:

```
conda activate assembly
```

### `--use-conda`

we are using rule-specific conda environments, because prokka and quast are annoying to install!

## wildcard constraints

this says: "the wildcard `lines` needs to be a number." c.f. regular expressions, and the [snakemake wildcard docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards).

this prevents it from trying to create `sample.full.out` with `head`.

### copy output contigs from assembly to rational name

```
# copy the final.contigs.fa file out of the megahit assembly directory
rule copy_genome_contigs:
    input:
        "{sample}_out"
    output:
        "{sample}.full.contigs.fa"
    shell:
        "cp {input}/final.contigs.fa {output}"
        
# copy the subset final.contigs.fa file out of the megahit assembly directory
rule copy_genome_contigs2:
    input:
        "{sample}.{lines}.out"
    output:
        "{sample}.{lines}.contigs.fa"
    shell:
        "cp {input}/final.contigs.fa {output}"
```

### run prokka to annotate the assembly

```
# annotate an assembly using prokka; *.contigs.fa -> *.prokka
rule annotate_contigs:
    input:
        "{prefix}.contigs.fa"
    output:
        directory("{prefix}.annot")
    threads: 8
    conda: "env-prokka.yml"
    shell:
        # note: a bug in prokka+megahit means we have to force success.
        # that's what "|| :" does.
        "prokka --outdir {output} --prefix {wildcards.prefix} {input} --cpus {threads} || :"
```

### run quast to evaluate the assembly

```

# evaluate an assembly using quast; *.contigs.fa => *.quast
rule quast_eval:
    input:
        "{prefix}.contigs.fa"
    conda: "env-quast.yml"      # note, you need to run this with --use-conda!
    output:
        directory("{prefix}.quast")
    shell: "quast {input} -o {output}"
```

## output of running

see the `.annot` directory for output of prokka, look at `.faa` for one example.

## k-mer abundance histograms

[sourmash](https://sourmash.readthedocs.io/en/latest/) is a tool produced by my lab for k-mer exploration.

Activate smash environment:
```
conda activate smash
```

```
sourmash sketch dna SRR2584857_?.fastq.gz --name SRR2584857 -o SRR2584857.full.sig -p abund
```

now run 

```
./kmer-abund-hist.py SRR2584857.full.sig --bins 20 --max 1000 --min 0
```

and add `--min` and `--max` until you get something like a bell curve :)

what is this doing?
