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
        self.tts_position = None  # None if no TTS has been found in a given region
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
        
        
def specify_gene_regions(gffutils_database : gffutils.FeatureDB, args : argparse.Namespace, name : str,
                         bacteria : str, samples : Dict[str, cf.ConfigSamplePredict]):
    """Specify regions for each gene in which we look for the TTS signatures. Returns TTS object with None value in place of the TTS position.
    Gene region is defined by args.tts-upstream, args.tts-downstream and args.tts-invasion.
    Gene regions shorter than 50nt are automatically excluded.

    Args:
        gffutils_database (gffutils.FeatureDB): gffutils database
        args (argparse.Namespace): argparse object
        name (str): sample id 
        bacteria (str): species name
        samples (Dict[str, cf.ConfigSamplePredict]): samples list
    """

    regions = []
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

        t = TTS(bacteria, [y for x, y in samples.items()], feature.seqid, feature.id, feature.strand, None, region_start,
                region_end, None)
        regions.append(t)
    return regions


def run_region_extraction(samples : Dict[str, cf.ConfigSample], gff_databases : dict, args : argparse.Namespace) -> List[TTS]:
    """Runs identification of the gene regions (in the form of TTS objects with no TTS)

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config
        gff_databases (Dict[str, gffutils.FeatureDB]:): dict storing references to gffutils dfatabases
        args (argparse.Namespace): argparse object

    Returns:
        list: list of TTS objects
    """
    
    regions = []  # list of tts objects (no tts found of course)
    for sample_id, sample_desc in samples.items():
        species = sample_desc.species
        annotations = sample_desc.annotations
        name = sample_desc.sample_id
        regions += specify_gene_regions(gff_databases[f"{annotations}.db"], args, name, species, samples)
    return regions