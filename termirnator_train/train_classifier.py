from argparse import ArgumentParser as Agp, FileType
from subprocess import Popen, PIPE
from os.path import join
from utils.configuration_file import ConfigSample
from utils.termination_sites import TTS, run_tts_extraction
from utils.parse_args import parse_arguments, check_species_validity
from utils.exceptions import SpeciesNotFoundInConfigFile
from utils.miscellaneous import eprint, downstream_gene
from utils.read_coverage import run_bam_to_bigWig_conversion
from utils.gene_annotations import create_gffutils_databases
from utils.create_datasets import create_training_set, create_testing_set, prepare_datasets
from utils.classifier import run_classifier_training
from Bio.Seq import Seq
from Bio import SeqIO
from typing import List, Dict, Tuple, Set
import numpy as np
import pandas as pd
import gffutils
import pyBigWig
import sys
import pybedtools
import pickle
import time
import os
import datetime
import argparse


def main():
    
    start_time = time.time()
    print("Reading config")
    args = parse_arguments()
    samples = ConfigSample.parse_config(args.config)
    
    # check if species provided in the command line are valid
    species_from_config = {x.species for x in samples.values()}
    species = args.test_on.strip().split(",")
    check_species_validity(species, species_from_config)
    species = args.train_on.strip().split(",")
    check_species_validity(species, species_from_config)
    
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"Converting BAMs to bigWigs. Running on {args.threads} threads")
    run_bam_to_bigWig_conversion(samples, args.threads)
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"Parsing annotation files")
    gff_databases = create_gffutils_databases(samples)
    chromosome_dict = make_chromosome_dict(samples)
    
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"TTS identification")
    tts = run_tts_extraction(samples, gff_databases, args)
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"Creating datasets. Running on {args.threads} threads")
    train_dataset, genes_to_keep_train = create_training_set(samples, tts, args, chromosome_dict, desired_size=args.sampling_size)
    test_dataset, genes_to_keep_test = create_testing_set(samples, tts, args, chromosome_dict, None)
    X_train, X_labels, Y_test, Y_labels, Y_info = prepare_datasets(train_dataset, test_dataset)

    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")

    start_time = time.time()
    print(f"Training and testing the classifier. Running on {args.threads} threads")
    classifier = run_classifier_training(X_train, X_labels, Y_test, Y_labels, Y_info, genes_to_keep_test, out_prefix=args.out_prefix,
                                         min_probability=args.min_prob,threads=args.threads)

    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    pybedtools.cleanup(remove_all=True)
    
    
if __name__ == "__main__":
    
    main()