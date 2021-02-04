---
tags: ggg, ggg2021, ggg298
---
# GGG 201(b), Lab Homework 1 - 2021

Due by 7pm, Tuesday Feb 9th

Please note that you may freely work together with others, but you must hand it in separately.

## 1. Sign up for GitHub Classroom

(Estimate: 10 min)

Sign up for GitHub Classroom using https://classroom.github.com/a/GJhfeDcK. This may involve creating a (free) GitHub account.

## 2. Log into farm & clone your github repo for hw1

(Estimate: 5 minutes)

(I will send you each an account name and password for farm, separately.)

On farm, clone your assignment repository, which will be something like `https://github.com/ngs-docs/2021-ggg-201b-lab-hw1-ctb`; do so like this,

```
git clone YOUR_REPO_URL 201b-lab-hw1
```

Change into that directory:
```
cd 201b-lab-hw1/
```
**You'll need to do this every time you log in.**

## 3. Create and activate the `vc` conda environment.

(Estimate: 15 minutes)

First, create the 'vc' conda environment, which is a collection of the
software you'll need in order to run the workflow:

```
conda env create -n vc -f environment.yml
```
**You only need to do this once.**

Then, activate the environment --
```
conda activate vc
```
**You'll need to run this every time you log in.**

and make sure you can successfully run the current workflow:

```
snakemake -j 1 --use-conda -p
```
which should produce a file `SRR2584857_1.ecoli-rel606.vcf`.

## 4. Edit and update your Snakefile

You'll need to edit your Snakefile in a text editor. If you're familiar with emacs or vi or some other UNIX command line editor, go for it! Just be sure that you're using spaces, and not tabs, when you indent lines. You can also use a local editor on your computer that supports remote SSH access, such as [Visual Studio Code](https://code.visualstudio.com/docs/remote/ssh).

If you don't have any prior experience here, I suggest using 'nano' -- see [this help document](https://www.redhat.com/sysadmin/getting-started-nano). I'll also put a video up shortly, demonstrating nano.

---

Adjust your Snakefile in two ways:

(1) as I showed in class, adjust your Snakefile to have a *separate* "unzip_genome" rule that creates an `ecoli-rel606.fa` file from `ecoli-rel606.fa.gz` file using `gunzip -c ecoli-rel606.fa.gz > ecoli-rel606.fa`.

(2) split the rule `call_variants` into three rules, one to call mpileup, one to run `bcftools call`, and one to run `bcftools view`.

Please make sure of the following:

* Running `snakemake -j 1 --use-conda -p` by itself should go from raw data to final results for all samples;
* All of the generated files (including intermediates) are in at least one 'output:' annotation, so that e.g. `snakemake --delete-all-output` removes all of the generated files in the directory.

## 5. Commit and push your changes back to github.

At any time, do the following to save changes.

```
git commit -am "updated Snakefile"
git push origin main
```
(you can do this as many times as you want, and save intermediate changes, etc. etc.)

One important note is that the only thing you should put in github are
the updated configuration files (Snakefile, in this case). All of the
generated files can be reproduced from the repo and so you don't need
to add them.

## 6. Relax in knowledge of a job well done.

(Please do inspect the Snakefile on github to make sure it has all your changes!)

## Questions? Problems?

If you run into any trouble with submission, that's ok - just let me know.

Please ask questions in the slack channel (be sure to tag in @ctitusbrown so I get notified!) or in the issue tracker for your personal github repository, which should https://github.com/ngs-docs/2021-ggg-201b-lab-hw1-USERNAME/issues (please tag in @ctb in the issue - this is Titus's github username).
