from argparse import ArgumentParser as Agp, FileType
from subprocess import Popen, PIPE
from os.path import join
from utils.configuration_file import ConfigSamplePredict
from utils.gene_regions_predict import TTS, run_region_extraction
from utils.parse_args import parse_arguments_predict
from utils.miscellaneous import make_chromosome_dict
from utils.read_coverage import run_bam_to_bigWig_conversion
from utils.gene_annotations import create_gffutils_databases
from utils.create_datasets_predict import create_datasets_predict, prepare_datasets
from utils.classifier import predict_using_classiffer
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
    args = parse_arguments_predict()
    samples = ConfigSamplePredict.parse_config(args.config)
    
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"Converting BAMs to bigWigs. Running on {args.threads} threads")
    run_bam_to_bigWig_conversion(samples, args.threads, which = ["rnaseq"])
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"Parsing annotation files")
    gff_databases = create_gffutils_databases(samples)
    chromosome_dict = make_chromosome_dict(samples)
    
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"TTS identification")
    tts = run_region_extraction(samples, gff_databases, args)
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    
    start_time = time.time()
    print(f"Creating datasets. Running on {args.threads} threads")
    dataset, genes_to_keep = create_datasets_predict(samples, tts, args, chromosome_dict)
    print(f"Finalizing datasets.")
    dataset, info = prepare_datasets(dataset)

    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")

    start_time = time.time()
    
    with open(args.classifier, "rb") as f:
        classifier = pickle.load(f)
    print(f"Training and testing the classifier. Running on {args.threads} threads")
    
    predict_using_classiffer(classifier, info, dataset, genes_to_keep, args.out_prefix)
 
    print(f"--- Execution time: {datetime.timedelta(seconds=time.time() - start_time)} ---")
    pybedtools.cleanup(remove_all=True)
    
    
if __name__ == "__main__":
    
    main()