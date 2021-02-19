---
tags: ggg, ggg2021, ggg201b
---
# Lab 5 outline - assembly! GGG 201(b) lab, Feb 5, 2021

## Can't log into farm?

You can run today's lab in binder, I think - 

[![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/binder-examples/r-conda/master?urlpath=rstudio)

and you will need to download these two files IF YOU ARE NOT ON FARM:

```
curl -L https://osf.io/aksmc/download -o SRR2584857_1.fastq.gz
curl -L https://osf.io/63pfr/download -o SRR2584857_2.fastq.gz
```

## Homework stuff, and farm stuff

homework updates:
- not online this weekend, sorry! I'll answer questions on monday morning
- push etc
- video, and farm updates (from annc)

### Running assembly on a farm node 

[description of farm](https://www.hpc.ucdavis.edu/farm-cluster) - we'll be using Farm III today, 'high2'.

Most HPCs or clusters at most universities distribute jobs across their cluster with a queuing system that allocations heterogeneous resources based on demand. This is a fancy way of saying that the queueing system tries to figure out how to efficiently use the available computers.

To do this, the cluster needs to know what resources you need. We'll use srun today to provide this, below.

On farm the queueing system we use is Slurm, and the way you access nodes is with 'srun'.

I should say that this is a quite different way to run programs from the way we've been doing it. A few differences:

* you run on a different computer than the head node, which means you have access to different (usually bigger!) resources (memory and disk and CPUs)
* the computer you run on may not have access to the Internet, though, so commands like `conda install` may not work properly.
* you run in a new environment, so you need to be sure to set your software up right (luckily conda makes this easy!)

In exchange for all of this complexity, you can set up dozens or hundreds or even thousands of jobs in advance, and then just ...walk away while slurm deals with it! (We'll talk about this in 298.)

### Some examples

- connect to farm
- activate conda
- log out/log back in, reactivate conda
- (am trying to avoid all of us connecting all at once ;)
- maybe do srun?

## Assembly! Roughly how it works...


reference dependent vs reference independent analyses
- reference based - has a genome => use genome to do stuff (alignment!)
- assembly based - don't have a genome => build a genome
- alignment free - doesn't need a genome => do analyses (typically k-mers)

(drawing / whiteboarding :)

### Running our first assembly!

NOTE: if you're on farm, schedule yourself a compute node rather than running on the head node --

```
srun --nodes=1 -p high2 -t 2:00:00 -c 4 --mem 6GB --pty /bin/bash
```
this says, "give me:
* one computer (`--nodes=1`)
* at high priority (`-p high2`)
* for 4 hours (`t 2:00:00`)
* with 4 CPUs (`-c 4`)
* and 6 GB of RAM (`--mem 6GB`)
* and set me up with an interactive login `--pty /bin/bash`)

Now, create yourself a subdirectory -

```
mkdir ~/lab5
cd ~/lab5
```

Next, we need to install the megahit software! I've already installed the [conda packaging system](https://github.com/ngs-docs/2021-GGG298/tree/latest/Week3-conda_for_software_installation) in your accounts, but we need to tell conda to use bioconda by adding optional *software channels* to it --

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

OK, now create a conda environment containing megahit:
```
conda create --name assembly megahit
```
and activate it:
```
conda activate assembly
```

Last, but not least, copy some data over:
```
cp ~ctbrown/data/ggg201b/SRR2584857_*.fastq.gz .
```
(on binder, use the curl at the top)

These are two files, approximately 190 MB each --
```
ls -lh
```

And... run!

```shell
megahit -1 SRR2584857_1.fastq.gz -2 SRR2584857_2.fastq.gz -f -m 5e9 -t 4
```

(Note here the `-m` is to limit memory usage to 5 GB, and the `-t` is to limit CPU usage. Otherwise megahit grabs ALL available memory and ALL available CPUs, which is kind of unfriendly in a shared environment!)

What does this output?

### Building a snakefile

We first need to install snakemake - do,

```
conda install -y snakemake-minimal
```

Now, let's build a Snakefile!

You can use `nano Snakefile` to edit the Snakefile locally.

To start that, let's make a single rule:
```
rule assemble:
    shell: """
        megahit -1 SRR2584857_1.fastq.gz -2 SRR2584857_2.fastq.gz -f -m 5e9 -t 4
    """
```
and you can run this with `snakemake -j 1 -p`.

### Log out!

If you're logged into farm, please type `exit` a few times to log out of your compute node - thanks!
