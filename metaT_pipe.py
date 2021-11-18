#!/usr/bin/env python3
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>
Date   : 2021-04-19
Purpose: Generate jobs scripts to be submitted to the UA HPC clusters
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Run modules for metatranscriptomics analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparser = parser.add_subparsers(help='Available modules')

    parser_cr = subparser.add_parser('create_reference',
                                     help='Module to create a reference for mapping')
    parser_cr.add_argument('input',
                           help='Input data to create reference')

    parser_cr.add_argument('--input_type',
                           help='Input data type: bins or contigs',
                           type=str,
                           choices=['contigs', 'bins'],
                           default='bins')

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Generate the job file"""

    args = get_args()


# --------------------------------------------------
if __name__ == '__main__':
    main()
