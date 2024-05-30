from typing import List, Dict, Tuple, Set
from . import configuration_file as cf
from subprocess import Popen, PIPE



def termseq_bam_to_bigWig(bam : str, threads : int):
    """Convert BAM files (mapped Term-seq reads) to bigWigs representing
    genomic coverage by 5' read ends ONLY
    (in Term-seq they represent transcript 3' termini).
    BigWig contains positions reflecting the first nt of the read.

    BamCoverage from deeptools is run to convert BAM to BigWig. Key options:
    --binSize 1 - 1-bp resolution
    --Offset 1 - use only the first base of each read to build the signal

    Args:
        bam (str): BAM file - mapped Term-seq reads
        threads (int): number of threads to use while running bamCoverage
        corresponding to genomic coverage by 5' read ends. Defaults to False.
    """
    p = Popen(["bamCoverage", "--bam", f"{bam}", "-o", f"{bam}.rev.bigwig",
               "--outFileFormat", "bigwig", "--binSize", "1", "--Offset",
               "1", "--normalizeUsing", "CPM", "--exactScaling",
               "--filterRNAstrand", "forward", '-p', f"{threads}"],
              stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    p = Popen(["bamCoverage", "--bam", f"{bam}", "-o", f"{bam}.fwd.bigwig",
               "--outFileFormat", "bigwig", "--binSize", "1", "--Offset",
               "1",  "--normalizeUsing", "CPM", "--exactScaling",
               "--filterRNAstrand", "reverse", '-p', f"{threads}"],
              stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
        

def rnaseq_bam_to_bigWig(bam : str, threads : int, library : str = "R"):
    """Convert BAM files (mapped RNA-seq reads) to bigWigs

    BamCoverage from deeptools is run to convert BAM to BigWig. Key options:
    --binSize 1 - 1-bp resolution
    --samFlagExclude 256 - exclude secondary alignments
    --normalizeUsing CPM - normalize coverage using CPM

    Args:
        bam (str): BAM file - mapped RNA-seq reads
        threads (int): number of threads to use while running bamCoverage
        library (str, optional): either R (reads corresponding to the gene on + strand will map to the - strand, and vice versa)
            or F (reads corresponding to the gene on + strand will map to the + strand, and vice versa)
    """
    if library == "R":
        st = "forward"
        st2 = "reverse"
    elif library == "F":
        st = "reverse"
        st2 = "forward"
    p = Popen(["bamCoverage", "--bam", f"{bam}", "-o", f"{bam}.rev.bigwig",
            "--outFileFormat", "bigwig", "--binSize", "1", "--normalizeUsing", "CPM", "--filterRNAstrand", st, '-p', f"{threads}"],
            stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    p = Popen(["bamCoverage", "--bam", f"{bam}", "-o", f"{bam}.fwd.bigwig",
            "--outFileFormat", "bigwig", "--binSize", "1", "--normalizeUsing", "CPM", "--filterRNAstrand", st2, '-p', f"{threads}"],
            stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    
def run_bam_to_bigWig_conversion(samples : Dict[str, cf.ConfigSample], threads : int, which : List[str] = ["rnaseq", "termseq"]):
    """Runs bam to bigWig conversion for each sample
    provided in the samples dictionary

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config
        threads (int): number of threads to use
        which (List[str]): which datasets to convert. Defaults to ["rnaseq", "termseq"]
    """
    
    for x in which:
        if x not in ['rnaseq', 'termseq']:
            raise ValueError(f"Value {x} is not in ['rnaseq', 'termseq']")
        
    if len(which) == 0:
        raise ValueError("Which parameter cannot be empty")

    for sample_id, sample_desc in samples.items():
        if 'termseq' in which:
            termseq_bam_to_bigWig(sample_desc.bam_path, threads)
        if 'rnaseq' in which:
            rnaseq_bam_to_bigWig(sample_desc.rna_seq_path, threads)