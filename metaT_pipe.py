#!/usr/bin/env python3
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>, Viviana Freire <vfreirezapata@email.arizona.edu>
Date   : 2021-04-19
Purpose: Generate jobs scripts to be submitted to the UA HPC clusters
"""

import argparse
import subprocess
import os
import glob
import pandas as pd


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Run modules for metatranscriptomics analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparser = parser.add_subparsers()

    parser_cr = subparser.add_parser('create_reference',
                                     help='Module to create a reference for mapping')
    parser_cr.add_argument('input_directory',
                           type=str,
                           help='Directory with the bins or contigs in fasta format'
                                'to create the reference for mapping')

    parser_cr.add_argument('--input_type',
                           type=str,
                           choices=['contigs', 'bins'],
                           help='Input data type: bins or contigs',
                           default='bins')

    parser_cr.add_argument('-o',
                           '--outdir',
                           type=str,
                           help='Output directory',
                           default='metaT_pipe_out')

    parser_cr.add_argument('-t',
                           '--threads',
                           help='Number of threads',
                           default=10)

    parser_cr.set_defaults(func=create_reference)

    parser_an = subparser.add_parser('annotate_reference',
                                     help='Module to annotate reference')

    parser_an.add_argument('input_reference',
                           type=str,
                           help='Directory with bins or contigs fasta file to annotate')

    parser_an.add_argument('reference_type',
                           type=str,
                           choices=['contigs', 'bins'],
                           help='Input data type: bins or contigs',
                           default='bins')

    parser_an.add_argument('--no_checkm',
                           action='store_true',
                           help='Do not run CheckM',
                           default=True)

    parser_an.add_argument('--no_gtdbtk',
                           action='store_true',
                           help='Do not run GTDB-TK',
                           default=True)

    parser_an.add_argument('-o',
                           '--outdir',
                           type=str,
                           help='Output directory',
                           default='metaT_pipe_out')

    parser_an.add_argument('-t',
                           '--threads',
                           help='Number of threads',
                           default=10)

    parser_an.set_defaults(func=annotate_reference)

    parser_map = subparser.add_parser('map_reads',
                                      help='Module to map reads to defined reference')

    parser_map.add_argument('--mapping_reference',
                            type=str,
                            help='Reference in fasta format to map the reads',
                            default=None)

    parser_map.add_argument('-r1',
                            type=str,
                            help='Forward reads for mapping in fastq format')

    parser_map.add_argument('-r2',
                            type=str,
                            help='Reverse reads for mapping in fastq format')

    parser_map.add_argument('--interleaved',
                            type=str,
                            help='Interleaved reads for mapping in fastq format')

    parser_map.add_argument('--mapper',
                            type=str,
                            choices=['bwa-mem', 'bowtie2'],
                            help='Choose to use either bwa-mem or bowtie2 for mapping the reads',
                            default='bwa-mem')

    parser_map.add_argument('-o',
                            '--outdir',
                            type=str,
                            help='Output directory',
                            default='metaT_pipe_out')

    parser_map.add_argument('-t',
                            '--threads',
                            help='Number of threads',
                            default=10)

    parser_map.set_defaults(func=map_reads)

    parser_gc = subparser.add_parser('get_read_counts',
                                     help='Module to obtain read counts from the mapping files')

    parser_gc.add_argument('mapping_directory',
                           type=str,
                           help='Directory with the mapping files in bam format')

    parser_gc.add_argument('--gff',
                           type=str,
                           required=True,
                           help='gff file with the gene coordinates to extract counts')

    parser_gc.add_argument('-o',
                           '--outdir',
                           type=str,
                           help='Output directory',
                           default='metaT_pipe_out')

    parser_gc.add_argument('-t',
                           '--threads',
                           help='Number of threads',
                           default=30)

    parser_gc.set_defaults(func=get_read_counts)

    args = parser.parse_args()

    # Check that command line arguments are specified properly

    return args


# --------------------------------------------------
def run_commands(cmd):
    """Function to run commands in the command line"""

    p = subprocess.run(
        cmd, shell=False, check=True, stderr=subprocess.STDOUT
    )

    return p


# --------------------------------------------------
def create_reference(args):
    """Create reference file for mapping"""

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    inputdir = os.path.join(args.input_directory, '*.fna')
    outdir = os.path.join(args.outdir, 'create_reference')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.input_type == 'bins':
        cmd = ['dRep', 'dereplicate', outdir, '-g', inputdir]
        run_commands(cmd)
    if args.input_type == 'contigs':
        cmd = ['cat', inputdir, '>', 'concat_contigs.fasta']
        run_commands(cmd)
        cmd = ['cd-hit-est', '-i', 'concat_contigs.fasta', '-o', 'temp_contigs99.fasta',
               '-c', '0.99', '-n', '10', '-T', str(args.threads)]
        run_commands(cmd)
        cmd = ['seqkit', 'seq', '-m', '2000', 'temp_contigs99.fasta', '>',
               os.path.join(outdir, 'derep_contigs99.fasta')]
        run_commands(cmd)
        cmd = ['rm', 'concat_contigs.fasta', 'temp_contigs99.fasta']
        run_commands(cmd)


# --------------------------------------------------
def annotate_reference(args):
    """Annotate reference"""

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    outdir = os.path.join(args.outdir, 'annotate_reference')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.reference_type == 'bins':
        if not args.no_checkm:
            cmd = ['checkm', 'lineage_wf', '-t', str(args.threads), '-x', 'fna', args.input_reference,
                   os.path.join(outdir, 'checkm_results')]
            run_commands(cmd)
            cmd = ['checkm', 'qa', os.path.join(outdir, 'checkm_results', 'lineage.ms'),
                   os.path.join(outdir, 'checkm_results'), '--tab_table', '-f',
                   os.path.join(outdir, 'checkm_results', 'checkm_table.tsv')]
            run_commands(cmd)
        if not args.no_gtdbtk:
            cmd = ['gtdbtk', 'classify_wf', '--genome_dir', args.input_reference, '--out_dir',
                   os.path.join(outdir, 'gtdb-tk_results'), '--cpus', str(args.threads)]
            run_commands(cmd)

        if args.no_checkm is False and args.no_gtdbtk is False:
            cmd = ['DRAM.py', 'annotate', '-i', "'" + os.path.join(args.input_reference, '*.fna') + "'", '-o',
                   os.path.join(outdir, 'dram_results'), '--checkm_quality',
                   os.path.join(outdir, 'checkm_results', 'checkm_table.tsv'), '--gtdb_taxonomy',
                   os.path.join(outdir, 'gtdb-tk_results', 'classify', 'gtdbtk.bac120.summary.tsv'),
                   '--threads', str(args.threads)]
            run_commands(cmd)

        elif args.no_checkm is True and args.no_gtdbtk is False:
            cmd = ['DRAM.py', 'annotate', '-i', "'" + os.path.join(args.input_reference, '*.fna') + "'", '-o',
                   os.path.join(outdir, 'dram_results'), '--gtdb_taxonomy',
                   os.path.join(outdir, 'gtdb-tk_results', 'classify', 'gtdbtk.bac120.summary.tsv'),
                   '--threads', str(args.threads)]
            run_commands(cmd)

        elif args.no_checkm is False and args.no_gtdbtk is True:
            cmd = ['DRAM.py', 'annotate', '-i', "'" + os.path.join(args.input_reference, '*.fna') + "'", '-o',
                   os.path.join(outdir, 'dram_results'), '--checkm_quality',
                   os.path.join(outdir, 'checkm_results', 'checkm_table.tsv'),
                   '--threads', str(args.threads)]
            run_commands(cmd)

        if args.no_checkm is True and args.no_gtdbtk is True:
            cmd = ['DRAM.py', 'annotate', '-i', "'" + os.path.join(args.input_reference, '*.fna') + "'", '-o',
                   os.path.join(outdir, 'dram_results'),
                   '--threads', str(args.threads)]
            run_commands(cmd)

    if args.reference_type == 'contigs':
        cmd = ['DRAM.py', 'annotate', '-i', args.input_reference, '-o',
               os.path.join(outdir, 'dram_results'),
               '--threads', str(args.threads)]
        run_commands(cmd)


# --------------------------------------------------
def map_reads(args):
    """Map reads to reference"""

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    outdir = os.path.join(args.outdir, 'map_reads')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.interleaved:
        cmd = ['coverm', 'make', '-r', args.mapping_reference, '--interleaved', args.interleaved, '-p', args.mapper,
               '-o', outdir, '-t', str(args.threads)]
        run_commands(cmd)
        bam_file = os.path.basename(args.mapping_reference) + '.' + os.path.basename(args.interleaved) + '.bam'
    else:
        cmd = ['coverm', 'make', '-r', args.mapping_reference, '-1', args.r1, '-2', args.r2, '-p', args.mapper,
               '-o', outdir, '-t', str(args.threads)]
        run_commands(cmd)
        bam_file = os.path.join(outdir, args.mapping_reference + args.interleaved + '.bam')

    filtered_bam_file = 'filtered' + bam_file
    cmd = ['coverm', 'filter', '-b', os.path.join(outdir, bam_file), '-o', os.path.join(outdir, filtered_bam_file), '-t', str(args.threads)]
    run_commands(cmd)

    sorted_bam_file = 'sorted.' + filtered_bam_file
    cmd = ['samtools', 'sort', os.path.join(outdir, filtered_bam_file), '-o', os.path.join(outdir, sorted_bam_file)]
    run_commands(cmd)


# --------------------------------------------------
def get_read_counts(args):
    """Get number of reads that mapped to each gene"""

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    outdir = os.path.join(args.outdir, 'get_read_counts')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    bam_files = glob.glob(args.mapping_directory + '**.bam')

    for file in bam_files:
        cmd = ['samtools', 'index', '-b', file]
        run_commands(cmd)
        cmd = ['dirseq', '--bam', file, '--gff', args.gff, '--measure_type', 'count', '>',
               os.path.join(outdir, os.path.basename(file) + '.counts.tsv')]
        run_commands(cmd)

    counts_table = pd.read_csv(glob.glob(outdir + '**.tsv')[0])[['ID']]
    counts_files = glob.glob(outdir + '**.tsv')
    for file in counts_files:
        temp = pd.read_csv(file, sep='\t')
        temp['final_count'] = temp['forward_read_count'] - temp['reverse_read_count']
        temp[temp['final_count'] < 0] = 0
        temp = temp[['ID', 'final_count']]
        counts_table = counts_table.merge(temp, on='ID')

    counts_table.to_csv(os.path.join(outdir, 'final_counts_table.csv'))


# --------------------------------------------------
def main():
    """Run metaT pipeline"""

    args = get_args()
    args.func(args)


# --------------------------------------------------
if __name__ == '__main__':
    main()
