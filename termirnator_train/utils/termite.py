from subprocess import Popen, PIPE
import os
from . import miscellaneous as msc

def run_TERMITe(replicates : list):
    """Finds transcription termination sites with TERMITe

    Args:
        replicates (list): list of replicates
    """
    
    bam_file_paths = [x.bam_path for x in replicates]
    rnaseq_bam_file_paths = [x.rna_seq_path for x in replicates]

    reverses = [x.replace(".bam", ".bam.fwd.bigwig") for x in bam_file_paths]
    forwards = [x.replace(".bam", ".bam.rev.bigwig") for x in bam_file_paths]
    
    rnaseq_reverses = [x.replace(".bam", ".bam.fwd.bigwig") for x in rnaseq_bam_file_paths]
    rnaseq_forwards = [x.replace(".bam", ".bam.rev.bigwig") for x in rnaseq_bam_file_paths]
    
    p = Popen(["mkdir", "-p", "tts_peaks"], stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    
    sample_names = [os.path.basename(x.replace(".sorted.bam", "")) for x in bam_file_paths]
    names = "_".join(sample_names)
    names_dash = "-".join(sample_names)
    

    p = Popen((["termite", "stable_rna_ends", "--rna-3prime-ends"] + forwards + ["--rna-seq-coverage"] + rnaseq_forwards + ["--sample-names"]
               + sample_names + ["--out-dir", "tts_peaks", "--strand", "forward", "--genome", f"{replicates[0].ref_genome}",
               "--name-prefix", f"{names}"]), stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = p.communicate()
    msc.eprint(stderr.decode("utf-8"))
    msc.eprint(stdout.decode("utf-8"))
    
    p = Popen((["termite", "stable_rna_ends", "--rna-3prime-ends"] + reverses + ["--rna-seq-coverage"] + rnaseq_reverses + ["--sample-names"]
               + sample_names + ["--out-dir", "tts_peaks", "--strand", "reverse", "--genome", f"{replicates[0].ref_genome}",
               "--name-prefix", f"{names}"]), stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = p.communicate()
    msc.eprint(stderr.decode("utf-8"))
    msc.eprint(stdout.decode("utf-8"))
    
    joined_forwards = ",".join(forwards)
    joined_reverses = ",".join(reverses)

    p = Popen((["termite", "annotate", "--termite-out-forward", f"tts_peaks/3prime_RNA_termini_{names_dash}_forward.narrowPeak",
                "--termite-out-reverse", f"tts_peaks/3prime_RNA_termini_{names_dash}_reverse.narrowPeak", "--sample-names",
                names, "--rna-3prime-ends-forward", f"{names}:{joined_forwards}",
                "--rna-3prime-ends-reverse", f"{names}:{joined_reverses}", "--trans-term-hp", "--rnafold",
                "--gene-annotations", f"{replicates[0].annotations}", "--genome", f"{replicates[0].ref_genome}",
                "--output", f"tts_peaks/{names}"]
               ), stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = p.communicate()
    msc.eprint(stderr.decode("utf-8"))
    msc.eprint(stdout.decode("utf-8"))

    
    
def extract_TTS_from_TERMITe_results(names : str):
    """extracts TTS from the TERMITe results and saves them in BED

    Args:
        names (str): replicate names joined with "_"
    """
    with open(f"tts_peaks/{names}.tsv") as termite_out:
        with open(f"final_peaks/{names}.bed","w") as out:
            for line in termite_out:
                sl = line.strip().split("\t")
                if sl[-1] == "+" and sl[-2] == "+":  # select only those with + in rnafold and transterhp column (intrinsic terminators)  
                    if "POT" not in line and "termite_id" not in line and "strand" not in line and "chromosome" not in line:
                        print(f"{sl[0]}\t{sl[3]}\t{sl[3]}\t{sl[6]}\t{sl[8]}\t{sl[5]}", file=out)