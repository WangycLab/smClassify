This script implements a complete preprocessing and classification workflow
for single-cell or bulk microbiome FASTQ data using Kraken2 + Bracken.


"""

import os
import re
import logging
import argparse
import subprocess as sp
from collections import defaultdict
from pyfastx import Fastq


def dir_check(dirname):
    """
    Create directory if it does not already exist.

    Parameters
    ----------
    dirname : str
        Path to directory to check / create.
    """
    if os.path.exists(dirname):
        return
    os.makedirs(dirname)


def file_check(filename):
    """
    Check whether a file exists.

    Parameters
    ----------
    filename : str
        Path to file.

    Returns
    -------
    bool
        True if file exists and is a regular file, else False.
    """
    return os.path.isfile(filename)


def read_config(config_file):
    """
    Read a simple key=value configuration file into a dictionary.

    Lines starting with '#' are treated as comments and ignored.
    Empty / whitespace-only lines are also ignored.

    Parameters
    ----------
    config_file : str
        Path to configuration file.

    Returns
    -------
    dict
        Dictionary mapping config keys to string values.
    """
    config = dict()
    with open(config_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            # Skip blank or whitespace-only lines
            if re.match(r'^\s$', line):
                continue
            # Strip whitespace, newline, etc.
            line = re.sub(r'[\r\n\s]+', '', line)
            arr = line.split('=')
            config[arr[0]] = arr[1]
    logging.info("Success: configuration info loaded.")
    return config


def run_kraken(fq_valid, config, suffix=''):
    """
    Run Kraken2 on the processed FASTQ file.

    Parameters
    ----------
    fq_valid : str
        Path to the final preprocessed FASTQ (or FASTQ.GZ) for classification.
    config : dict
        Configuration dictionary, must contain keys:
        'sample', 'outdir', 'kraken', 'krakenDb', and optionally 'process'.
    suffix : str, optional
        Optional suffix to append to output file names (e.g. for batch runs).

    Returns
    -------
    str
        Path to the Kraken2 report file.
    """
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
    """
    Run Bracken on a Kraken2 report to refine taxonomic abundance estimates.

    Parameters
    ----------
    config : dict
        Configuration dictionary, must contain keys:
        'sample', 'outdir', 'krakenDb', 'braken', and optionally
        'process', 'filter_threshold'.
    kreport : str
        Path to the Kraken2 report file (*.report).
    suffix : str, optional
        Optional suffix to append to output file names.

    Notes
    -----
    - Bracken is run separately at species level (-l S) and genus level (-l G).
    - Classified reads from Kraken2 are passed to Bracken using '-f'.
    """
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
    # Pass classified FASTQ to Bracken if present (plain or gzipped)
    if file_check(fq_classify):
        cmds += ['-f', fq_classify]
    elif file_check(fq_classify + ".gz"):
        # For gzipped input, use process substitution with zcat
        cmds += ['-f', "<(zcat " + fq_classify + ".gz)"]
    else:
        print("Error: missing classified fastq file.")
        exit(1)

    # Species level
    cmd_S = cmds + ['-l', 'S', '-o', output]
    # Genus level
    cmd_G = cmds + ['-l', 'G', '-o', output_g]
    sp.run(" ".join(cmd_S), shell=True, stdout=fh)
    sp.run(" ".join(cmd_G), shell=True, stdout=fh)

    fh.close()


class PreFlight:
    """
    PreFlight: handle all preprocessing steps before running Kraken/Bracken.

    Responsibilities:
        - Set up project directory structure (clean_data, Result, etc.).
        - Combine paired-end FASTQ reads and track per-barcode read counts.
        - Select top-N cell barcodes based on read depth.
        - Deduplicate reads using an external deduplication program.
        - Retrieve final FASTQ for selected barcodes using helper binaries.

    Methods
    -------
    run():
        Execute the full preflight pipeline and return the path to final FASTQ.
    pre_dedup():
        Combine FASTQ reads and compute cell barcode read counts.
    dedup(prededup):
        Run external deduplication tool on combined FASTQ.
    ready_fq(fq_dedup, cellNum):
        Extract reads for selected barcodes and generate the final FASTQ.
    """

    def __init__(self, config):
        """
        Initialize a PreFlight instance from configuration.

        Parameters
        ----------
        config : dict
            Configuration dictionary including keys:
            'sample', 'outdir', 'nubeam_dedup', and optionally 'bclen', 'mode'.
        """
        self.config = config
        self.sample = self.config['sample']
        proj_dir = self.config['outdir']

        # Data directory to hold intermediate and cleaned FASTQ files
        self.data_dir = os.path.join(proj_dir, 'clean_data')
        dir_check(self.data_dir)

        # Result directory for output files
        self.res_dir = os.path.join(proj_dir, 'Result')
        dir_check(self.res_dir)

    def run(self):
        """
        Run the full preflight pipeline:
            1) Combine R1/R2 and build barcode read counts
            2) Deduplicate reads
            3) Retrieve reads for top-N barcodes

        Returns
        -------
        str
            Path to the final preprocessed FASTQ file.
        """
        prededup, cellNum = self.pre_dedup()
        fq_dedup = self.dedup(prededup)
        fq_valid = self.ready_fq(fq_dedup, cellNum)
        return fq_valid

    def pre_dedup(self):
        """
        Pre-deduplication step: combine R1 and R2, track barcode counts,
        and select the threshold for top-N barcodes (default N = 200,000).

        Returns
        -------
        tuple
            (prededup_fastq_path, cellNum)

            prededup_fastq_path : str
                Path to the combined pre-deduplicated FASTQ file.
            cellNum : int
                Number of barcodes above the selected read-count threshold.
        """
        read_counts = defaultdict(int)
        fq1 = os.path.join(self.data_dir, '{}_1.fq'.format(self.sample))
        fq2 = os.path.join(self.data_dir, '{}_2.fq'.format(self.sample))

        # Prefer plain FASTQ; fall back to gzipped if necessary
        fq1 = fq1 if file_check(fq1) else fq1 + ".gz"
        fq2 = fq2 if file_check(fq2) else fq2 + ".gz"
        if not file_check(fq1) or not file_check(fq2):
            print("Missing fastq files.")
            exit(1)

        bc_len = self.config.get('bclen', 20)
        prededup = os.path.join(self.data_dir, '{}.prededup.fq'.format(self.sample))

        # Build pre-dedup FASTQ and barcode read counts if not already present
        if not file_check(prededup):
            fh = open(prededup, 'w')
            for item1, item2 in zip(
                Fastq(fq1, build_index=False, full_name=True),
                Fastq(fq2, build_index=False, full_name=True)
            ):
                id_1, seq_1, qual_1 = item1
                id_2, seq_2, qual_2 = item2

                # Filter out very short R2 reads (likely low-quality)
                if len(seq_2) < 60:
                    continue

                cb = seq_1[0:bc_len]
                read_counts[cb] += 1
                new_seq = seq_1 + seq_2
                new_qual = qual_1 + qual_2

                fh.write('@{}\n{}\n+\n{}\n'.format(id_2, new_seq, new_qual))
            fh.close()

        x, y = [], []
        n = 1
        reads_ctx = os.path.join(self.data_dir, '{}_read_counts.tsv'.format(self.sample))

        # If read count file exists, re-use it instead of recomputing
        if file_check(reads_ctx):
            read_counts = defaultdict(int)
            with open(reads_ctx, 'r') as infile:
                for line in infile:
                    arr = line.strip().split('\t')
                    read_counts[arr[0]] = int(arr[1])

        # Sort barcodes by read counts (descending)
        sorted_read_counts = sorted(read_counts.items(), key=lambda item: item[1], reverse=True)

        # Determine read-count threshold for top 200,000 barcodes (or fewer if not enough)
        num_thres = min(200000, len(sorted_read_counts)) - 1
        countthres = sorted_read_counts[num_thres][1]

        cellNum = 0
        if file_check(reads_ctx):
            # If counts file is present, just count how many barcodes are above threshold
            for cb, val in sorted_read_counts:
                if val >= countthres:
                    cellNum += 1
        else:
            # Otherwise, write counts to file and compute cellNum
            with open(reads_ctx, 'w') as fh:
                for cb, val in sorted_read_counts:
                    if val >= countthres:
                        cellNum += 1
                    fh.write('{}\t{}\n'.format(cb, val))
                    y.append(val)
                    x.append(n)
                    n += 1
        return prededup, cellNum

    def dedup(self, prededup):
        """
        Run external deduplication on the pre-deduplicated FASTQ file.

        Parameters
        ----------
        prededup : str
            Path to the pre-deduplicated FASTQ file.

        Returns
        -------
        str
            Path to the deduplicated FASTQ file.
        """
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

    def ready_fq(self, fq_dedup, cellNum):
        """
        Prepare final FASTQ for Kraken2 by retrieving reads from selected barcodes.

        Parameters
        ----------
        fq_dedup : str
            Path to deduplicated FASTQ file.
        cellNum : int
            Number of barcodes to keep (top-N based on read counts).

        Returns
        -------
        str
            Path to the final FASTQ file with selected cells.
        """
        fq_final = os.path.join(self.data_dir, '{}_final_2.fq'.format(self.sample))
        reads_ctx = os.path.join(self.data_dir, '{}_read_counts.tsv'.format(self.sample))
        prefix = os.path.join(self.data_dir, self.sample)
        dir_name = os.path.abspath(os.path.dirname(__file__))
        mode = self.config.get('mode', 'm20')

        # Select correct helper binary for barcode-based read retrieval
        if mode == 'm20':
            reads_retriev = os.path.join(dir_name, 'readsRetriev2_tmp')
        elif mode == '10x':
            reads_retriev = os.path.join(dir_name, 'readsRetriev3')
        elif mode == 'test':
            reads_retriev = os.path.join(dir_name, 'readsRetriev2_test')
        else:
            # If an unknown mode is provided, still try to use 'readsRetriev2_tmp'
            reads_retriev = os.path.join(dir_name, 'readsRetriev2_tmp')

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
    """
    Entry point for the microbiome classification pipeline.

    Command-line options
    --------------------
    --cfg <FILE>           : Path to configuration file (required).
    -k / --kraken-only     : Only run Kraken2 (skip Bracken).
    -b / --bracken-only    : Only run Bracken on an existing Kraken2 report.
    -s / --suffix <STR>    : Optional suffix for output file names.
    """
    AP = argparse.ArgumentParser(
        description="microbiome classification pipeline.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    AP.add_argument('--cfg', metavar='', help='The config file', required=True)
    AP.add_argument('-k', '--kraken-only', dest='kraken', action='store_true',
                    help="Only run kraken")
    AP.add_argument('-b', '--bracken-only', dest='bracken', action='store_true',
                    help="Only run braken")
    AP.add_argument('-s', '--suffix', type=str, default='',
                    help="Suffix for kraken output files. (Default: '')")

    args = AP.parse_args()
    config = read_config(args.cfg)

    # Bracken-only mode: expect an existing Kraken2 report
    if args.bracken:
        kreport = os.path.join(
            config['outdir'],
            'Result',
            '{}_kraken{}.report'.format(config['sample'], args.suffix)
        )
        if file_check(kreport):
            run_braken(config, kreport, args.suffix)
            return
        else:
            print("Cannot find kraken report file. Need to rerun from the beginning.")

    preflight = PreFlight(config)

    fq_valid = os.path.join(preflight.data_dir, '{}_final_2.fq'.format(preflight.sample))

    # If final FASTQ already exists (plain or gzipped), reuse it;
    # otherwise, rerun the preflight pipeline.
    if file_check(fq_valid + ".gz"):
        fq_valid = fq_valid + ".gz"
    elif not file_check(fq_valid):
        fq_valid = preflight.run()

    # Run Kraken2 classification
    kreport = run_kraken(fq_valid, config, args.suffix)

    # If user requested Kraken-only mode, do not run Bracken
    if args.kraken:
        return

    # Run Bracken post-classification refinement
    run_braken(config, kreport, args.suffix)


if __name__ == '__main__':
    # Basic logging config (INFO level is sufficient for pipeline messages)
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s'
    )
    main()
