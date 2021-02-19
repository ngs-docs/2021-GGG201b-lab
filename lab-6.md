---
tags: ggg, ggg2021, ggg201b
---
# Lab 6 outline - assembly, part 2! GGG 201(b) lab, Feb 12, 2021

### Preparing for today if you used binder last week

**you only need to do this if you didn't do it last week**

That having been said, you can run the commands below just fine :)

Log into farm, and then do:

```
mkdir ~/lab5
cd ~/lab5

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create --name assembly megahit snakemake-minimal
```

Then, copy some data over:

```
cp ~ctbrown/data/ggg201b/SRR2584857_*.fastq.gz ./
```

### Revisiting: running assembly on a farm node

Log into farm.

Schedule yourself a compute node rather than running on the head node --

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

Now, go to your existing subdirectory - you may need to go up to "preparing for today" if you get an error that `lab5` doesn't exist.

```
cd ~/lab5
```

Finally, activate your `assembly` conda environment:
```
conda activate assembly
```

### A first assembly Snakefile

Last week, we created a small Snakefile with a single rule:
```
rule assemble:
    shell: """
        megahit -1 SRR2584857_1.fastq.gz -2 SRR2584857_2.fastq.gz -f -m 5e9 -t 4
    """
```
which can be run with

```
snakemake -p -j 1
```

The problem is this takes a few minutes to run, and everytime we run it it's going to run again! Let's fix up the Snakefile a bit -- add an `output:` line, like so:

```
    output: directory("megahit_out")
```

Here, the `directory(...)` is telling snakemake that this directory is _completely_ under snakemake's control, and that it can safely delete it & recreate everything under it with running that rule.

Let's also add some inputs --

```
    input:
        "SRR2584857_1.fastq.gz",
        "SRR2584857_2.fastq.gz",
```

### Adding wildcards and changing the directory name

(Will do interactively, and paste final Snakefile here)

### Add a default rule at the top

at the top, add

```
rule all:
    input: "SRR2584857.out"
```

### RESTING POINT:

```
rule assemble_all_samples:
    input:
        "SRR2584857_out"

rule assemble:
    input:
        r1 = "{sample}_1.fastq.gz",
        r2 = "{sample}_2.fastq.gz",
    output: directory("{sample}_out")
    shell: """
        megahit -1 {input.r1} -2 {input.r2} -f -m 5e9 -t 4 -o {output}
    """
```

### Calculating genome coverage

To estimate coverage, we need to know

(1) how many reads
(2) how long each read is (on average)
(3) how big the genome is

and then multiply reads by their length, and divide by the size of the genome.

So - how big is the genome, approximately? How would we find out? (I'll show you a bad way to estimate it, computationally :)

To count the reads:
```
gunzip -c SRR2584857_1.fastq.gz | wc -l
```
\- note there have to be exactly that many reads in the \_2 file as well.

How do you figure out how long each read is?

Try:

```
gunzip -c SRR2584857_1.fastq.gz | head
```

(Note FastQC will tell you all of this read information, too.)

### Subsampling for coverage exploration

Let's try subsampling the reads a bit!

Write a subsampling rule:

```
rule subsample_reads:
    input:
        r1="{sample}_1.fastq.gz",
        r2="{sample}_2.fastq.gz",
    output:
        r1="{sample}_1.{lines}.fastq.gz",
        r2="{sample}_2.{lines}.fastq.gz",
    shell: """
        gunzip -c {input.r1} | \
                head -{wildcards.lines} | \
                gzip > {output.r1} \
                || true
        gunzip -c {input.r2} | \
                head -{wildcards.lines} | \
                gzip > {output.r2} \
                || true
    """
```

(There's ...a lot of magic in here, sorry!)

To explain a few things in the shell block -
* the \\ at the end of the lines is a *continuation* marker that says "this is really all one command". Importantly, it cannot be followed by anything except a return (newline).
* so, there are only two command lines in there.
* they are basically the same, so let's just look at one...
* the \| is a *pipe*, which says, "feed the output of the command before the pipe into the input of the command after the pipe"
* in order, the pipes say,
    * uncompress the read file and send to output
    * take the first `{wildcards.lines}` number of lines - where is it getting this from? - and send to output
    * compress the result and send to file `{output.r1}`
    * the || says, if the previous commands fail, it's ok - this is something snakemake needs! (demonstrate, titus!)
    * (this is needed because the gunzip fails when 'head' says, ENOUGH)

So: this rule takes some number of lines and puts it in a file.

* why are there two wild cards?
* how do we run this rule?! Ask snakemake for `SRR2584857_1.100000.fastq.gz`!

### How do we modify the assembly rule to assemble the subsets?

Well, first, it is simplest to copy it and modify it...

```
rule assemble_subset:
    input:
        r1="{sample}_1.{lines}.fastq.gz",
        r2="{sample}_2.{lines}.fastq.gz",
    output: directory("{sample}.{lines}.out")
    shell: """
        megahit -1 {input.r1} -2 {input.r2} -f -m 5e9 -t 4 -o {output}
    """
```

and then let's ask snakemake to run it for us!

```
snakemake -j 1 -p SRR2584857.100000.out
```

### Add to top of Snakefile

In-class task: make your input line for the default rule at the top of the Snakefile assemble subsets 100000 and 200000 in size.

### RESTING POINT

```
rule DEBUG_assemble_all_samples:
    input:
        "SRR2584857.100000.out",

rule GO_WALKABOUT:
    input:
        "SRR2584857.200000.out",
        "SRR2584857.300000.out",
        "SRR2584857.500000.out",

rule CRUSH_MY_ENEMIES:
    input: "SRR2584857.10000000.out"

rule assemble:
    input:
        r1 = "{sample}_1.fastq.gz",
        r2 = "{sample}_2.fastq.gz",
    output: directory("{sample}_out")
    shell: """
        megahit -1 {input.r1} -2 {input.r2} -f -m 5e9 -t 4 -o {output}
    """

rule assemble_subset:
    input:
        r1 = "{sample}_1.{lines}.fastq.gz",
        r2 = "{sample}_2.{lines}.fastq.gz",
    output: directory("{sample}.{lines}.out")
    shell: """
        megahit -1 {input.r1} -2 {input.r2} -f -m 5e9 -t 4 -o {output}
    """
    
# some random comment
rule subsample_reads:
    message: "Subsampling {wildcards.sample} to {wildcards.lines} lines"
    input:
        r1="{sample}_1.fastq.gz",
        r2="{sample}_2.fastq.gz",
    output:
        r1="{sample}_1.{lines}.fastq.gz",
        r2="{sample}_2.{lines}.fastq.gz",
    shell: """
        # do this thing for r1
        gunzip -c {input.r1} | \
                head -{wildcards.lines} | \
                gzip > {output.r1} \
                || true
        gunzip -c {input.r2} | \
                head -{wildcards.lines} | \
                gzip > {output.r2} \
                || true
    """
```

### A few notes --

This is "parameter exploration" in that we can explore many different subset sizes quite easily!

Note also that megahit runs WAY faster for smaller subsets...

...so this is a pretty useful way to test snakefiles!

(It's what I did with the variant calling pipeline, actually :)

## End of lab

### Log out!

If you're logged into farm, please type `exit` a few times to log out of your compute node - thanks!
