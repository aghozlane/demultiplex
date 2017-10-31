#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
from __future__ import print_function
import os
import sys
import argparse
#import pandas as pd
import gzip
import csv

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2017, Institut Pasteur"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-i1', dest='fastq_file_R1', type=isfile,
                        help='Fastq R1 file to demultiplex.')
    parser.add_argument('-i2', dest='fastq_file_R2', type=isfile,
                        help='Fastq R2 file to demultiplex.')
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        help='Single end Fastq file to demultiplex.')
    parser.add_argument('-a', dest='assignation_tab_file', type=isfile,
                        required=True, help='Table assigning index to sample.')
    #parser.add_argument('-t', dest='tag_fastq', type=str, default='',
    #                    help='Tag in filename for paired reads.')
    parser.add_argument('-o', dest='output_dir', default= os.curdir + os.sep,
                        type=isdir, help='Output directory (default ./).')
    args = parser.parse_args()
    return args


def load_assign(assignation_tab_file):
    """
    """
    run_dict = {}
    try:
        with open(assignation_tab_file, "rt") as assign:
            assign_reader = csv.reader(assign, delimiter="\t")
            for line in assign_reader:
                run_dict[line[0]] = line[1]
    except IOError:
        sys.exit("Error cannot open {0}".format(assignation_tab_file))
    return run_dict


def demultiplex(fastq_file, run_dict, output_dir, tag=""):
    """
    """
    assigned_count = 0
    notassigned_count = 0
    total_count = 0
    output_dict = {}
    assigned_dict = {}
    notassigned_list = []
    try:
        #notassigned_filename = (output_dir + os.sep + "not_assigned"
        #                        + tag_fastq + ".fastq.gz")
        #notassigned_output = gzip.open(notassigned_filename, "wb")
        ext = os.path.splitext(fastq_file)[1]
        if ext == ".gz":
            fastq = gzip.open(fastq_file, "rb")
        else:
            fastq = open(fastq_file, "rt")
        for index in run_dict.keys():
            output_filename = (output_dir + os.sep + run_dict[index]
                                   + tag + ".fastq.gz")
            output_dict[index] = gzip.open(output_filename, "wb")
        # Start reading
        for line in fastq:
            #print(line)
            header = line
            sequence = fastq.next()
            tag = fastq.next()
            quality = fastq.next()
            #print(sequence)
            #print(sequence[0:8])
            #key = run_df['Index']['Sample'] == sequence[0:8]
            #for
            total_count += 1
            if sequence[0:8] in run_dict:
                output = output_dict[sequence[0:8]]
                #print(output_filename)
                output.write("{0}{1}{2}{3}".format(
                    header, sequence[8:], tag, quality[8:]))
                assigned_count += 1
                assigned_dict[header.split(' ')[0]] = sequence[0:8]
            else:
                ##TODO More accurate search
                #print(sequence)
                #notassigned_output.write("{0}{1}{2}{3}".format(header, sequence,
                #                                               tag, quality))
                notassigned_count += 1
                notassigned_list.append(header.split(' ')[1])
                #print("Not found: {0}".format(sequence[0:8]))
                #for index in run_dict.keys():
                #    if index in sequence:
                #        print("found {0}".format(index))

        # Close file
        fastq.close()
        #notassigned_output.close()
        for index in run_dict.keys():
            output_dict[index].close()
        print("Reads assigned: {0} not assigned: {1} total: {2}".format(
                assigned_count, notassigned_count, total_count))
    except IOError:
        sys.exit("Error cannot open {0}".format(fastq_file))
    return assigned_dict, notassigned_list


def assign_R2(fastq_file, assigned_dict, run_dict, output_dir):
    """
    """
    output_dict = {}
    try:
        ext = os.path.splitext(fastq_file)[1]
        if ext == ".gz":
            fastq = gzip.open(fastq_file, "rb")
        else:
            fastq = open(fastq_file, "rt")
        for index in run_dict.keys():
            output_filename = (output_dir + os.sep + run_dict[index]
                                   + '_R2' + ".fastq.gz")
            output_dict[index] = gzip.open(output_filename, "wb")
        for line in fastq:
            header = line
            sequence = fastq.next()
            tag = fastq.next()
            quality = fastq.next()
            if header.split(' ')[0] in assigned_dict:
                output = output_dict[assigned_dict[header.split(' ')[0]]]
                output.write("{0}{1}{2}{3}".format(
                    header, sequence[8:], tag, quality[8:]))
        # Close file
        fastq.close()
        for index in run_dict.keys():
            output_dict[index].close()
    except IOError:
        sys.exit("Error cannot open {0}".format(fastq_file))


def main():
    """Main program
    """
    args = getArguments()
    # Load assignation
    #run_df = pd.read_csv(args.assignation_tab_file, sep="\t", names=['index', 'sample'])
    run_dict = load_assign(args.assignation_tab_file)
    print("Keys:")
    print(run_dict)
    if args.fastq_file_R1 and args.fastq_file_R2:
        # Demultiplex
        assigned_dict, notassigned_list = demultiplex(
            args.fastq_file_R1, run_dict, args.output_dir, '_R1')
        # Assign R1 based on R2
        assign_R2(args.fastq_file_R2, assigned_dict, run_dict, args.output_dir)
    elif args.fastq_file:
        demultiplex(args.fastq_file, run_dict, args.output_dir)
    else:
        sys.exit("Error demultiplexing is impossible to perform : "
                 "no file provided or one pair is missing")


if __name__ == '__main__':
    main()
