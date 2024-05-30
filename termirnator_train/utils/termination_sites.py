import os
import argparse
from subprocess import Popen, PIPE
import gffutils
from . import termite
from . import miscellaneous
from . import configuration_file as cf
from typing import List, Dict, Tuple, Set


class TTS():
    def __init__(self, bacteria, samples, chromosome, gene, strand, tts_position, region_start, region_end,
                 average_height):
        self.bacteria = bacteria
        self.samples = samples
        self.chromosome = chromosome
        self.gene = gene
        self.tts_position = tts_position  # can be None if no TTS has been found in a given region
        self.strand = strand
        self.region_start = region_start
        self.region_end = region_end
        self.average_height = average_height
        
    def __eq__(self, x):
        if self.bacteria == x.bacteria and self.samples == x.samples and self.chromosome == x.chromosome \
            and self.gene == x.gene and self.tts_position == x.tts_position and self.strand == x.strand:
            return True
        else:
            return False
        
        
def specify_gene_regions(gffutils_database : gffutils.FeatureDB, args : argparse.Namespace, names : str):
    """Specify regions for each gene in which we look for the TTS signatures. Saves regions in BED.
    Gene region is defined by args.tts-upstream, args.tts-downstream and args.tts-invasion.
    Gene regions shorter than 50nt are automatically excluded.

    Args:
        gffutils_database (gffutils.FeatureDB): gffutils database
        args (argparse.Namespace): argparse object
        names (str): replicate names joined with "_"
    """
    with open(f"tts_peaks/gene_regions_{names}.bed", "w") as f:
        for feature in gffutils_database.features_of_type("gene"):

            # specify region of the gene
            # in which to look for TTS signatures

            if feature.strand == "+":
                closest_downstream = miscellaneous.downstream_gene(gffutils_database,
                                                    feature.seqid, "+",
                                                    feature.end)
                region_start = feature.end - args.tts_upstream
                region_end = feature.end + args.tts_downstream
                region_end = (min(region_end, closest_downstream+args.tts_invasion)
                            if closest_downstream is not None else region_end)
                region_end = (feature.end + args.tts_invasion
                            if region_start >= region_end else region_end)
            elif feature.strand == "-":
                closest_downstream = miscellaneous.downstream_gene(gffutils_database,
                                                    feature.seqid, "-",
                                                    feature.start)
                region_start = feature.start - args.tts_downstream
                region_end = feature.start + args.tts_upstream
                region_start = (max(region_start - args.tts_invasion, closest_downstream)
                                if closest_downstream is not None
                                else region_start)
                region_start = (feature.start - args.tts_invasion 
                                if region_start >= region_end else region_start)

            if abs(region_end - region_start) < 50:
                continue

            f.write(f"{feature.seqid}\t{max(0, region_start)}\t{region_end}\t{feature.id}\t0\t{feature.strand}\n")


def extract_replicates(sample : cf.ConfigSample,
                       parsed_config : Dict[str, cf.ConfigSample]) -> List[str]:
    """Extract all replicates of a given sample
    provided in parse_config results

    Args:
        sample (ConfigSample): sample
        parsed_config (Dict[str, ConfigSample]): dictinary returned by parse_config

    Returns:
        List[ConfigSample]: list of ConfigSample objects
    """
    
    res_samples = []

    replicate = parsed_config[sample.sample_id].condition
    for sample_id, sample_desc in parsed_config.items():
        if sample_desc.condition == replicate:
            res_samples.append(sample_desc)
    return res_samples


def extract_TTS_from_gene_regions(bacteria : str, names : str, replicates : list) -> List[TTS]:
    """Extracts TTS found within the specified gene regions (see specify gene regions).
    Selects only the most abundant TTS found within the gene region

    Args:
        bacteria (str): species name
        names (str): replicate names joined with "_"
        replicates (list): list of replicates
        
    Returns:
        list: list of TTS objects
    """

    tts = []
    
    with open(f"tts_peaks/{names}.intersect.bed","w") as out:
            p = Popen((["bedtools", "intersect", "-a", f"tts_peaks/gene_regions_{names}.bed", "-b", f"final_peaks/{names}.bed",
                        "-wao", "-s"]),
                    stdout=out, stderr=PIPE)
            stdout, stderr = p.communicate()
            miscellaneous.eprint(stderr.decode("utf-8"))
        

    with open(f"tts_peaks/{names}.intersect.bed") as f:
        intersect = {}  ## used to select only the most abundant termination site for each gene
        for line in f:
            sl = line.strip().split("\t")
            gene = sl[3]
            prominence = float(sl[10])
            if gene in intersect and prominence > intersect[gene][0]:  
                intersect[gene] = (prominence, sl)
            elif gene not in intersect:
                intersect[gene] = (prominence, sl)

        for gene, value in intersect.items():
            sl = value[1]
            if sl[7] == "-1":
                termination_site = TTS(bacteria, replicates, sl[0],
                                sl[3], sl[5], None, int(sl[1]), int(sl[2]), sl[10])
                tts.append(termination_site)
            else:
                termination_site = TTS(bacteria, replicates, sl[0],
                                sl[3], sl[5], int(sl[7]), int(sl[1]), int(sl[2]), sl[10])
                tts.append(termination_site)
    miscellaneous.eprint(f"Number of TTS for {bacteria} ({replicates}): {len([x for x in tts if x.tts_position is not None])} ")
    
    return tts
        
        
def extract_replicates_for_TTS(tts: TTS, samples : Dict[str, cf.ConfigSample]) -> List[cf.ConfigSample]:
    """Returns samples that correspond to the TTS

    Args:
        tts (TTS): TTS object
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config

    Returns:
        List[ConfigSample]: list of samples
    """
    
    result = []
    for sample, sample_desc in samples.items():
        if sample_desc.bam_path in [x.bam_path for x in tts.samples]:
            result.append(sample_desc)
    return result


def extract_tts(bacteria : str, replicates : list, gffutils_database : gffutils.FeatureDB, path : str,
                args : argparse.Namespace) -> List[TTS]:
    """Function extracts TTS positions from regions downstream of the genes
    annotated in the gffutils_database

    Args:
        bacteria (str): species name
        replicates (list): list of replicates
        gffutils_database (gffutils.FeatureDB): Connection to the database
        path (str): path to bam file
        args (argparse.Namespace): argparse object

    Returns:
        list: list of TTS objects
    """
    
    termite.run_TERMITe(replicates)

    p = Popen(["mkdir", "-p", "final_peaks"], stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    
    p = Popen(["mkdir", "-p", "tts_peaks"], stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    
    bam_file_paths = [x.bam_path for x in replicates]
    rnaseq_bam_file_paths = [x.rna_seq_path for x in replicates]
    sample_names = [os.path.basename(x.replace(".sorted.bam", "")) for x in bam_file_paths]
    names = "_".join(sample_names)
    
    specify_gene_regions(gffutils_database, args, names)
    
    termite.extract_TTS_from_TERMITe_results(names)
    
    tts = extract_TTS_from_gene_regions(bacteria, names, replicates)

    return tts


def run_tts_extraction(samples : Dict[str, cf.ConfigSample], gff_databases : dict, args : argparse.Namespace) -> List[TTS]:
    """Runs TTS identification

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config
        gff_databases (Dict[str, gffutils.FeatureDB]:): dict storing references to gffutils dfatabases
        args (argparse.Namespace): argparse object

    Returns:
        list: list of TTS objects
    """
    
    done_samples = []
    tts = []  # list of identified termination events
    for sample_id, sample_desc in samples.items():
        if sample_desc.bam_path not in done_samples:
            replicates = extract_replicates(sample_desc, samples)
            done_samples += [x.bam_path for x in replicates]
            species = sample_desc.species
            annotations = sample_desc.annotations
            path = sample_desc.bam_path
            tts += extract_tts(species, replicates, gff_databases[f"{annotations}.db"], path, args)
    return tts