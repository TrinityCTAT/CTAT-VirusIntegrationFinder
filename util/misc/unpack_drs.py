#!/usr/bin/env python3

import sys, os, re
import glob
import subprocess
import logging


def main():

    usage = "\n\n\tusage: {} sample_id fastqs.tar[.gz]\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    sample_id = sys.argv[1]
    fastqs_tar = sys.argv[2]

    os.makedirs("fastq")
    
    subprocess.check_call(f"tar -xvf {fastqs_tar} -C fastq", shell=True)

    fq_files = sorted(glob.glob("fastq/*"))

    def is_gzipped(filename):
        if re.search(".gz$", filename):
            return True
        else:
            return False
    
    if len(fq_files) == 1:
        if is_gzipped(fq_files[0]):
            os.rename(fq_files[0], sample_id + "_1.fq.gz")
        else:
            os.rename(fq_files[0], sample_id + "_1.fq")
            subprocess.check_call("gzip " + sample_id + "_1.fq", shell=True)

    elif len(fq_files) == 2:
        if is_gzipped(fq_files[0]):
            os.rename(fq_files[0], sample_id + "_1.fq.gz")
            os.rename(fq_files[1], sample_id + "_2.fq.gz")
        else:
            os.rename(fq_files[0], sample_id + "_1.fq")
            os.rename(fq_files[1], sample_id + "_2.fq")
            subprocess.check_call(f"gzip {sample_id}_*.fq", shell=True)

    elif len(fq_files) == 4:
        method = "cat"
        if is_gzipped(fq_files[0]):
            method = "zcat"

        subprocess.check_call(f"{method} {fq_files[0]} {fq_files[1]} | gzip -c > {sample_id}_1.fq.gz", shell=True)
        subprocess.check_call(f"{method} {fq_files[2]} {fq_files[3]} | gzip -c > {sample_id}_2.fq.gz", shell=True)

        


if __name__=='__main__':
    main()

    
