import argparse
import multiprocessing
from typing import Dict, Tuple, List
import pandas as pd
from . import configuration_file as cf
from . import termination_sites as ts
import pyBigWig
import pybedtools
from Bio.Seq import Seq
import random
import concurrent.futures
from statistics import median

import warnings
warnings.filterwarnings("ignore")



def filter_by_avg_exp(dataset : pd.DataFrame, return_all_genes : bool = False) -> Dict[str, int]:
    """Function that filters genes having the average RNA-seq coverage of their first 20 nucleotides above median.
    Returns the dict of gene identifiers meeting this criteria.

    Args:
        dataset (pd.DataFrame): dataset
        return_all_genes (bool, optional): return all genes. Do not filter. Defaults to False.

    Returns:
        Dict[str, int]: dictionary storing gene names to keep (keys) and their average expression levels (first 20nt, values) 
    """
    
    genes_to_keep = {}
    
    if return_all_genes:
        for index, row in dataset.iterrows():
            genes_to_keep[row["gene_x_sample"]] = 0
        return genes_to_keep

    data_temp = dataset.groupby('gene_x_sample', as_index=False).first()
    data = data_temp.loc[: , "cov_0":"cov_20"]
    data['average'] = data.mean(axis=1)
    means = list(data['average'])
    data["gene_x_sample"] = data_temp["gene_x_sample"]
    data = data.loc[data['average'] >= median(means)]
    
    for index, row in data.iterrows():
        genes_to_keep[row["gene_x_sample"]] = row["average"]
    del data
    
    return genes_to_keep


def create_dataset(samples : Dict[str, cf.ConfigSample], tts_list : List[ts.TTS], upstream : int,
                   downstream : int, species : List[str], chromosome_dict : Dict[str, str],
                   genes_with_no_tts : bool = False,
                   ) -> pd.DataFrame:
    """Creates training or testing dataset for machine learning models

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config
        tts_list (List[TTS]): list of TTS objects corresponding to the identified 3' RNA termini
        upstream (int): extract up to \{u\} nucleotides upstream from the 3\' RNA termini
        downstream (int): extract up to \{u\} nucleotides downstream from the 3\' RNA termini
        species (List[str]): list of species (same as in the config file) to be included in the dataset
        genes_with_no_tts (bool): whether to keep the gene regions without TTS

    Returns:
        pd.DataFrame: training/testing dataset
    """
    
    d = {}
    for tts in tts_list:
        for sample in tts.samples:  #! replaced the code above
            if sample.species not in species:
                continue
            
            if tts.tts_position == None and genes_with_no_tts == False:
                continue  # skip genes with no clear TTS
            if tts.strand == "+":
                bw = pyBigWig.open(f"{sample.rna_seq_path}.rev.bigwig")
            elif tts.strand == "-":
                bw = pyBigWig.open(f"{sample.rna_seq_path}.fwd.bigwig")
                
            seq_size = bw.chroms(tts.chromosome)
            coverage = bw.values(tts.chromosome, 0, seq_size)
            seq = chromosome_dict[sample.ref_genome][tts.chromosome]
            
            region_len = tts.region_end - tts.region_start
            for position in range(region_len):
                pos = position + tts.region_start
                
                if tts.strand == "+":
                    start, end = pos - upstream, pos + downstream
                else:
                    start, end = pos - downstream, pos + upstream

                sequence = seq[start-1:end]
                sequence = sequence if tts.strand == "+" else str(Seq(sequence).reverse_complement())

                cov = coverage[start-1:end]
                cov = cov[::-1] if tts.strand == "-" else cov
                
                if len(cov) != abs(end-start)+1 or abs(end-start)+1 != len(sequence):
                    break
                
                if "gene" in d.keys():
                    d["gene"].append(tts.gene)
                else:
                    d["gene"] = [tts.gene]
                
                if "label" in d.keys():
                    d["label"].append(pos == tts.tts_position)
                else:
                    d["label"] = [pos == tts.tts_position]

                if "genomic_position_chr" in d.keys():
                    d["genomic_position_chr"].append(tts.chromosome)
                else:
                    d["genomic_position_chr"] = [tts.chromosome]

                if "genomic_position_pos" in d.keys():
                    d["genomic_position_pos"].append(pos)
                else:
                    d["genomic_position_pos"] = [pos]

                if "genomic_position_strand" in d.keys():
                    d["genomic_position_strand"].append(tts.strand)
                else:
                    d["genomic_position_strand"] = [tts.strand]

                if "rnaseq_sample" in d.keys():
                    d["rnaseq_sample"].append(sample.rna_seq_path)
                else:
                    d["rnaseq_sample"] = [sample.rna_seq_path]

                for i, c in enumerate(cov):
                    if f"cov_{i}" in d.keys():
                        d[f"cov_{i}"].append(c)
                    else:
                        d[f"cov_{i}"] = [c]

                for i, nuc in enumerate(sequence):
                    if f"seq_{i}" in d.keys():
                        d[f"seq_{i}"].append(nuc)
                    else:
                        d[f"seq_{i}"] = [nuc]
    try:
        df = pd.DataFrame(d)
    except Exception as e:
        print(e)
    
    return df


def filter_and_subsample_dataset(df : pd.DataFrame, desired_size : int = None) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Runs filtering by the expression level and sampling

    Args:
        df (pd.DataFrame): dataset
        desired_size (int, optional): if specified, limit the number of positive and negative examples to this value by sampling. Defaults to None.

    Returns:
        Tuple[pd.DataFrame, Dict[str, float]]: dataset, genes to be kept
    """

    df["gene_x_sample"] = df["gene"] + df["rnaseq_sample"]
    genes_to_keep = filter_by_avg_exp(df, return_all_genes=False)
    df = df.loc[df['gene_x_sample'].isin(genes_to_keep.keys())]
    del df["gene_x_sample"]
    
    # sampling
    if desired_size:
        positive_examples = df.loc[df['label'] == True]
        negative_examples = df.loc[df['label'] == False]
        try:
            negative_examples_sampled = negative_examples.sample(desired_size)
        except ValueError:
            negative_examples_sampled = negative_examples.sample(desired_size, replace=True)
        try:
            positive_examples_sampled = positive_examples.sample(desired_size)
        except ValueError:
            positive_examples_sampled = positive_examples.sample(desired_size, replace=True)
        
        df = pd.concat([negative_examples_sampled, positive_examples_sampled], axis=0)
        del negative_examples_sampled
        del positive_examples_sampled
    
    return df, genes_to_keep
    


def prepare_datasets(train_dataset : pd.DataFrame, test_dataset : pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame,
                                                                                                       pd.DataFrame, pd.DataFrame,
                                                                                                       pd.DataFrame]:
    """Prepare training and testing datasets for RandomForestClassifier (subsampling, filtering genes with insufficient expression)

    Args:
        train_dataset (pd.DataFrame): training dataset
        test_dataset (pd.DataFrame): testing dataset

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]: training examples,
            training labels, testing examples, testing labels, testing dataset info
    """
    
    seq_cols = [x for x in list(train_dataset.columns) if ("seq_" in x and x.startswith("seq"))]
    train_dataset = pd.get_dummies(train_dataset, columns=seq_cols, drop_first=True)
    train_dataset = train_dataset.sample(frac=1).reset_index(drop=True)
    train_labels = train_dataset.iloc[:, 1].values
    train_dataset = train_dataset.iloc[:, 6:].values
    
    seq_cols = [x for x in list(test_dataset.columns) if ("seq_" in x and x.startswith("seq"))]
    test_dataset = pd.get_dummies(test_dataset, columns=seq_cols, drop_first=True)
    test_dataset = test_dataset.reset_index(drop=True)
    test_labels = test_dataset.iloc[:, 1].values
    test_info = test_dataset.iloc[:, :7]
    test_dataset = test_dataset.iloc[:, 6:].values
    
    return train_dataset, train_labels, test_dataset, test_labels, test_info



def create_training_set(samples : Dict[str, cf.ConfigSample], tts_list : List[ts.TTS], args : argparse.Namespace,
                        chromosome_dict : Dict[str, str], desired_size : int = None) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Runs create_datasets

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config
        tts_list (List[TTS]): list of TTS objects corresponding to the identified 3' RNA termini
        args (argparse.Namespace): argparse object
        desired_size (int, optional): if specified, limit the number of positive and negative examples to this value by sampling. Defaults to None.

    Returns:
        Tuple[pd.DataFrame, Dict[str, float]]: training dataset, genes to be kept (with expression level of the examined 20nt above the median)
    """
    
    divided_tts_list = divide_tts_list(tts_list)
    
    species = args.train_on.strip().split(",")
    
    create_dataset_args = [(samples, x, args.upstream, args.downstream, species, False) for x in divided_tts_list.values()]
    
    args_samples = [samples for x in divided_tts_list.values()]
    args_x = [x for x in divided_tts_list.values()]
    args_upstream = [args.upstream for x in divided_tts_list.values()]
    args_downstream = [args.downstream for x in divided_tts_list.values()]
    args_species = [species for x in divided_tts_list.values()]
    args_false = [False for x in divided_tts_list.values()]
    args_chromosome_dicts = [chromosome_dict for x in divided_tts_list.values()]
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(create_dataset, args_samples, args_x, args_upstream,
                               args_downstream, args_species, args_chromosome_dicts, args_false)
        
    train_dataset = pd.concat(results)
        
    # filtering by expression levels and subsampling
    train_dataset, genes_to_keep = filter_and_subsample_dataset(train_dataset, desired_size)
    
    return train_dataset, genes_to_keep



def create_testing_set(samples : Dict[str, cf.ConfigSample], tts_list : List[ts.TTS], args : argparse.Namespace,
                       chromosome_dict : Dict[str, str], desired_size : int = None) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Runs create_datasets

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config
        tts_list (List[TTS]): list of TTS objects corresponding to the identified 3' RNA termini
        args (argparse.Namespace): argparse object
        desired_size (int, optional): if specified, limit the number of positive and negative examples to this value by sampling. Defaults to None.

    Returns:
        Tuple[pd.DataFrame, Dict[str, float]]: testing dataset, genes to be kept (with expression level of the examined 20nt above the median)
    """
    
    divided_tts_list = divide_tts_list(tts_list)
    
    species = args.test_on.strip().split(",")
    
    create_dataset_args = [(samples, x, args.upstream, args.downstream, species, True) for x in divided_tts_list.values()]
    
    
    args_samples = [samples for x in divided_tts_list.values()]
    args_x = [x for x in divided_tts_list.values()]
    args_upstream = [args.upstream for x in divided_tts_list.values()]
    args_downstream = [args.downstream for x in divided_tts_list.values()]
    args_species = [species for x in divided_tts_list.values()]
    args_false = [True for x in divided_tts_list.values()]
    args_chromosome_dicts = [chromosome_dict for x in divided_tts_list.values()]
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(create_dataset, args_samples, args_x, args_upstream,
                               args_downstream, args_species, args_chromosome_dicts, args_false)
        
    test_dataset = pd.concat(results)

    
    # filtering by expression levels and subsampling
    test_dataset, genes_to_keep = filter_and_subsample_dataset(test_dataset, desired_size)
    
    return test_dataset, genes_to_keep


def divide_tts_list(tts_list : List[ts.TTS]) -> Dict[Tuple[str], List[ts.TTS]]:
    """Divide TTS list into different lists that can be run in bulk on multiple threads

    Args:
        tts_list (List[TTS]): list of TTS objects

    Returns:
        Dict[Tuple[str], List[TTS]]: divided lists of TTS saved as a dict values
    """
    div_list = {}
    for tts in tts_list:
        samples = tuple(tts.samples)
        if samples not in div_list:
            div_list[samples] = [tts]
        else:
            div_list[samples].append(tts)
            
    return div_list



