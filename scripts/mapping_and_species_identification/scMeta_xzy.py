#!/usr/bin/env python

import os
import re
import logging
import argparse
import subprocess as sp
from collections import defaultdict
from pyfastx import Fastq

def dir_check(dirname):
    if os.path.exists(dirname):
        return
    os.makedirs(dirname)

def file_check(filename):
    return os.path.isfile(filename)

def read_config(config_file):
    config = dict()
    with open(config_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if re.match(r'^\s$', line):
                continue
            line = re.sub(r'[\r\n\s]+','',line)
            arr = line.split('=')
            config[arr[0]] = arr[1]
    logging.info("Success: configuration info loaded.")
    return config

def run_kraken(fq_valid, config, suffix=''):
    sample = config['sample']
    res_dir = os.path.join(config['outdir'], 'Result')
    kraken2 = config['kraken']
    kdb = config['krakenDb']

    fq_classify = os.path.join(res_dir, 'output_classified{}.fq'.format(suffix))
    fq_unclassify = os.path.join(res_dir, 'output_unclassified{}.fq'.format(suffix))

    output = os.path.join(res_dir, '{}_kraken{}.output'.format(sample, suffix))
    report = os.path.join(res_dir, '{}_kraken{}.report'.format(sample, suffix))

    sp.run([
        kraken2,
        '--threads', config.get('process', '8'),
        '--db', kdb,
        '--unclassified-out', fq_unclassify,
        '--classified-out', fq_classify,
        '--report', report,
        '--output', output,
        fq_valid
    ])    
    return report

def run_braken(config, kreport, suffix=''):
    sample = config['sample']
    res_dir = os.path.join(config['outdir'], 'Result')
    output = os.path.join(res_dir, '{}_sc_taxonomy{}.report'.format(sample, suffix))
    output_g = os.path.join(res_dir, '{}_sc_taxonomy_G{}.report'.format(sample, suffix))
    b_log = os.path.join(res_dir, '{}.braken{}.log'.format(sample, suffix))
    kmer_db = os.path.join(config['krakenDb'], 'database100mers.kmer_distrib')
    fq_classify = os.path.join(res_dir, 'output_classified{}.fq'.format(suffix))
    cutoff = config.get('filter_threshold', '0')

    cmds = [
        config['braken'],
        '-i', kreport,
        '-k', kmer_db,
        '-t', '0',
        '-p', config.get('process', '8'),
        '-c', cutoff
    ]

    fh = open(b_log, 'w')
    if file_check(fq_classify):
        cmds += ['-f', fq_classify]
    elif file_check(fq_classify + ".gz"):
        cmds += ['-f', "<(zcat " + fq_classify + ".gz)"]
    else:
        print("Error: missing classified fastq file.")
        exit(1)

    cmd_S = cmds + ['-l', 'S', '-o', output]
    cmd_G = cmds + ['-l', 'G', '-o', output_g]
    sp.run(" ".join(cmd_S), shell=True, stdout=fh)
    sp.run(" ".join(cmd_G), shell=True, stdout=fh)

    # sp.run([
    #     config['braken'],
    #     '-i', kreport,
    #     '-o', output,
    #     '-k', kmer_db,
    #     '-f', fq_classify,
    #     '-l', 'S',
    #     '-t', '0',
    #     '-p', config.get('process', '8'),
    #     '-c', cutoff
    # ], stdout=fh)

    # sp.run([
    #     config['braken'],
    #     '-i', kreport,
    #     '-o', output_g,
    #     '-k', kmer_db,
    #     '-f', fq_classify,
    #     '-l', 'G',
    #     '-t', '0',
    #     '-p', config.get('process', '8'),
    #     '-c', cutoff
    # ], stdout=fh)

    fh.close()


class PreFlight:
    def __init__(self, config):
        self.config = config
        self.sample = self.config['sample']
        proj_dir = self.config['outdir']

        self.data_dir = os.path.join(proj_dir, 'clean_data')
        dir_check(self.data_dir)

        self.res_dir = os.path.join(proj_dir, 'Result')
        dir_check(self.res_dir)  

    def run(self):
        prededup,cellNum = self.pre_dedup()
        fq_dedup = self.dedup(prededup)
        fq_valid = self.ready_fq(fq_dedup, cellNum)
        return fq_valid

    def pre_dedup(self):

        read_counts = defaultdict(int)
        fq1 = os.path.join(self.data_dir, '{}_1.fq'.format(self.sample))
        fq2 = os.path.join(self.data_dir, '{}_2.fq'.format(self.sample))
        
        # if not file_check(fq1) or not file_check(fq2):
        #     print("Missing fastq files.")
        #     exit(1)

        fq1 = fq1 if file_check(fq1) else fq1 + ".gz"
        fq2 = fq2 if file_check(fq2) else fq2 + ".gz"
        if not file_check(fq1) or not file_check(fq2):
            print("Missing fastq files.")
            exit(1)

        # mode = self.config.get('mode', 'm20')
        bc_len = self.config.get('bclen', 20)

        prededup = os.path.join(self.data_dir, '{}.prededup.fq'.format(self.sample))

        if not file_check(prededup):
            fh = open(prededup, 'w')
            for item1, item2 in zip(Fastq(fq1, build_index=False, full_name=True), \
                Fastq(fq2, build_index=False, full_name=True)):
                id_1, seq_1, qual_1 = item1
                id_2, seq_2, qual_2 = item2

                if len(seq_2) < 60:
                    continue

                # if mode == 'm20':
                #     cb = seq_1[0:20]
                # elif mode == '10x':
                #     cb = seq_1[0:16]
                # elif mode == 'test':
                #     cb = seq_1[0:21]

                cb = seq_1[0:bc_len]
                read_counts[cb] += 1
                new_seq = seq_1 + seq_2
                new_qual = qual_1 + qual_2

                fh.write('@{}\n{}\n+\n{}\n'.format(id_2, new_seq, new_qual))
            fh.close()

        x, y = [], []
        n = 1
        reads_ctx = os.path.join(self.data_dir, '{}_read_counts.tsv'.format(self.sample))

        if file_check(reads_ctx):
            read_counts = defaultdict(int)
            with open(reads_ctx, 'r') as infile:
                for line in infile:
                    arr = line.strip().split('\t')
                    read_counts[arr[0]] = int(arr[1])

        sorted_read_counts = sorted(read_counts.items(), key=lambda item: item[1], reverse=True)
        ### get the readcounts thershold (top 200000)
        # countthres = sorted_read_counts[199999][1]
        num_thres = min(200000, len(sorted_read_counts)) - 1
        countthres = sorted_read_counts[num_thres][1]
        cellNum = 0
        if file_check(reads_ctx):
            for cb, val in sorted_read_counts:
                if val >= countthres:
                    cellNum += 1
        else:
            with open(reads_ctx, 'w') as fh:
                for cb, val in sorted_read_counts:
            #   for cb, val in sorted(read_counts.items(), key=lambda item: item[1], reverse=True):
                    if val >= countthres:
                        cellNum += 1
                    fh.write('{}\t{}\n'.format(cb, val))
                    y.append(val)
                    x.append(n)
                    n += 1
        return prededup, cellNum

    def dedup(self, prededup):
        fq_dedup = os.path.join(self.data_dir, '{}.dedup.fq'.format(self.sample))
        if file_check(fq_dedup):
            return fq_dedup
        nubeam_dedeup = self.config['nubeam_dedup']

        sp.run([
            nubeam_dedeup,
            '-i', prededup,
            '-o', fq_dedup
        ])

        return fq_dedup
              
    #def ready_fq(self, fq_dedup, cell_num=200000):
    def ready_fq(self, fq_dedup,cellNum):
        fq_final = os.path.join(self.data_dir, '{}_final_2.fq'.format(self.sample))
        reads_ctx = os.path.join(self.data_dir, '{}_read_counts.tsv'.format(self.sample))
        prefix = os.path.join(self.data_dir, self.sample)
        dir_name = os.path.abspath(os.path.dirname(__file__))
        mode = self.config.get('mode', 'm20')
        # reads_retriev = os.path.join(dir_name, 'readsRetriev2_tmp')

        if mode == 'm20':
            reads_retriev = os.path.join(dir_name, 'readsRetriev2_tmp')
        elif mode == '10x':
            reads_retriev = os.path.join(dir_name, 'readsRetriev3')
        elif mode == 'test':
            reads_retriev = os.path.join(dir_name, 'readsRetriev2_test')

        sp.run([
            reads_retriev,
            '-bc', reads_ctx,
            '-cells', str(cellNum),
            '-fq', fq_dedup,
            '-o', prefix,
            '-tab', 'on'
            ])

        return fq_final

def main():
    AP = argparse.ArgumentParser(
            description="microbiome classification pipeline.",
            formatter_class=argparse.RawTextHelpFormatter,
        )
    
    AP.add_argument('--cfg', metavar='', help='The config file', required=True)
    AP.add_argument('-k', '--kraken-only', dest='kraken', action='store_true', help="Only run kraken")
    AP.add_argument('-b', '--bracken-only', dest='bracken', action='store_true', help="Only run braken")
    AP.add_argument('-s', '--suffix', type=str, default='', help="Suffix for kraken output files. (Default: '')")

    args = AP.parse_args()
    config = read_config(args.cfg)

    if args.bracken:
        kreport = os.path.join(config['outdir'], 'Result', '{}_kraken{}.report'.format(config['sample'], args.suffix))
        if file_check(kreport):
            run_braken(config, kreport, args.suffix)
            return
        else:
            print("Cannot find kraken report file. Need to rerun from the beginning.")

    preflight = PreFlight(config)

    fq_valid = os.path.join(preflight.data_dir, '{}_final_2.fq'.format(preflight.sample))

    ## Check if need to rerun dedup
    if file_check(fq_valid + ".gz"):
        fq_valid = fq_valid + ".gz"
    elif not file_check(fq_valid):
        fq_valid = preflight.run()

    # if not file_check(fq_valid):
    #     if file_check(fq_valid + ".gz"):
    #         fq_valid = fq_valid + ".gz"
    #     else:
    #         fq_valid = preflight.run()

    # if not file_check(fq_valid):
    #     fq_valid = preflight.run()

    kreport = run_kraken(fq_valid, config, args.suffix)
#    fq_valid = preflight.run()
    if args.kraken:
        return
    
    run_braken(config, kreport, args.suffix)

if __name__ == '__main__':
    main()

