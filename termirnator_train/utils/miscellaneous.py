import sys
import gffutils
from Bio import SeqIO


def eprint(*args, **kwargs):
    """Prints to stderr
    """
    print(*args, file=sys.stderr, **kwargs)


def downstream_gene(gffutils_database : gffutils.FeatureDB, chromosome : str,
                    strand : str, gene_end : int) -> int:
    """Get the nearest downstream gene 5' end coordinate from the gene_end.
    Must be on the same strand.

    Args:
        gffutils_database (gffutils.FeatureDB): Connection to the database
        chromosome (str): chromosome/contig ID
        strand (str): strand (+/-)
        gene_end (int): position from which we start looking for the nearest
                       downstream 5' gene end

    Returns:
        int: position of the nearest downstream 5' gene end
    """

    if strand == "+":
        cmd = 'SELECT min(start), id from features where seqid == ' \
            '"{}" AND strand == \"+\" AND ' \
            'featuretype == \"gene\" AND end > {}'.format(chromosome, gene_end)
        downstream = gffutils_database.execute(cmd)
    else:
        cmd = 'SELECT max(end), id from features where seqid == ' \
            '"{}" AND strand == \"-\" AND featuretype == \"gene\" ' \
            'AND start < {}'.format(chromosome, gene_end)
        downstream = gffutils_database.execute(cmd)
    for i in downstream.fetchone():
        return i
    
    
def make_chromosome_dict(samples):
    chromosome_dict = {}
    genomes = [x.ref_genome for x in samples.values()]
    for g in genomes:
        fasta_sequences = SeqIO.parse(open(g),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if g in chromosome_dict:
                chromosome_dict[g][name] = sequence
            else:
                chromosome_dict[g] = {name : sequence}
    return chromosome_dict
