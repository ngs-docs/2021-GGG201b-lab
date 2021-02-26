#! /usr/bin/env python
import sys
import argparse
import subprocess
import csv
import pprint


def main():
    p = argparse.ArgumentParser()
    p.add_argument('samples', nargs='+')
    p.add_argument('-S', '--sensitive', action='store_true')
    args = p.parse_args()

    store_type = "specific"
    if args.sensitive:
        store_type = "sensitive"

    for sample in args.samples:
        vcf = f"{sample}.{store_type}.vcf"
        bam = f"{sample}.bam.sorted"
        genome = "SRR2584857_1.ecoli-rel606.fa"

        variants = []
        with open(vcf, 'rt') as fp:
            r = csv.reader(fp, delimiter='\t')
            for row in r:
                if row[0][0] == '#':
                    continue
                variants.append(row)


        for variant in variants:
            acc = variant[0]
            position = variant[1]
            ref = variant[3]
            alt = variant[4]

            view_pos = int(position) - 20
            view_pos = max(0, view_pos)

            print(position, ref, alt)
            subprocess.check_call(f"samtools index {bam}", shell=True)
            cmd = f"samtools tview -p {acc}:{view_pos} {bam} {genome}"
            print(cmd)
            subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    sys.exit(main())
