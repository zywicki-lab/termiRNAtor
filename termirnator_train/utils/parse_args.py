# -*- coding: utf-8 -*-
from argparse import ArgumentParser as Agp
from argparse import FileType
from argparse import ArgumentTypeError
from . import exceptions

from typing import List, Set

def check_positive(value : int):
    ivalue = int(value)
    if ivalue <= 0:
        raise ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def parse_arguments():
    """Parse commandline arguments

    Returns:
        argparse.Namespace: argparse.Namespace object
    """
    parser = Agp(description='Train and test termiRNAtor classifier')
    parser.add_argument('-t', '--term-seq-samples', type=FileType('r'),
                        required=True, help='config file containing \
                            information about term-seq samples. Template can \
                            be found in config.csv', dest='config')
    parser.add_argument('-p', '--threads', type=check_positive,
                        required=False, default=1, help='Number of threads to \
                        use', dest='threads')
    parser.add_argument('-u', '--upstream', type=int, help='extract up to \{u\}\
                                                            nucleotides upstream\
                                                            from the tested position\
                                                            as a training/testing\
                                                            feature',
                        required=False, default=150, dest="upstream")                        
    parser.add_argument('-d', '--downstream', type=int, help='extract up to \{d\}\
                                                            nucleotides downstream\
                                                            from the tested position\
                                                            as a training/testing\
                                                            feature',
                        required=False, default=20, dest="downstream")
    parser.add_argument('--tts-downstream', type=int, help='search for the TTS up to {--tts-downstream}\
                                                            from the gene 3\' end',
                        required=False, default=150, dest="tts_downstream")
    parser.add_argument('--tts-upstream', type=int, help='search for the TTS up to {--tts-upstream}\
                                                            from the gene 3\' end',
                        required=False, default=20, dest="tts_upstream")
    parser.add_argument('--tts-invasion', type=int, help='search for the TTS in the region with maximal overlap of {--tts-invasion} nt\
                                                            with the nearest downstream gene',
                        required=False, default=10, dest="tts_invasion")
    parser.add_argument('--train_on', type=str, help='comma-separated list of organisms\
                                                      to train the model (same as\
                                                      species in the config file)',
                        required=False, default="", dest="train_on")
    parser.add_argument('--test_on', type=str, help='comma-separated list of organisms\
                                                      to test the model (same as\
                                                      species in the config file)',
                        required=True, dest="test_on")
    parser.add_argument('--training_set_size', type=int, default=2000000,
                        help='Desired sampling size for the training dataset',
                        required=False, dest="sampling_size")
    parser.add_argument('--min_probability', type=float, default=.15,
                        help='min probability of the classifier output to be considered a 3\'end',
                        required=False, dest="min_prob")
    parser.add_argument('--out-prefix', type=str, help='prefix of the output files',
                        required=False, default="termiRNAtor", dest="out_prefix")
    args = parser.parse_args() 

    return args


def check_species_validity(species : List[str], species_from_config : Set[str]):
    """Checks if species names proided with --train_on and 
    --test_on flags are valid and present in the config file

    Args:
        species (List[str]): list of the species provided with flags
        species_from_config (Set[str]): all species present in the config file

    Raises:
        exceptions.SpeciesNotFoundInConfigFile: at least one species provided in the\
            command line does not exist in the config file
    """
    for s in species:
        if s not in species_from_config:
            raise exceptions.SpeciesNotFoundInConfigFile(s)


def parse_arguments_predict():
    """Parse commandline arguments for the prediction module

    Returns:
        argparse.Namespace: argparse.Namespace object
    """
    parser = Agp(description='Predict TTS based on RNA-seq data')
    parser.add_argument('-r', '--rna-seq-samples', type=FileType('r'),
                        required=True, help='config file containing \
                            information about rna-seq samples. Template can \
                            be found in config_prediction.csv', dest='config')
    parser.add_argument('-c', '--classifier', type=str,
                        required=True, help='Pickled classifier being an output\
                            of the train_classifier module', dest='classifier')
    parser.add_argument('-p', '--threads', type=check_positive,
                        required=False, default=1, help='Number of threads to \
                        use', dest='threads')
    parser.add_argument('-u', '--upstream', type=int, help='extract up to \{u\}\
                                                            nucleotides upstream\
                                                            from the tested position\
                                                            as a training/testing\
                                                            feature',
                        required=False, default=150, dest="upstream")                        
    parser.add_argument('-d', '--downstream', type=int, help='extract up to \{d\}\
                                                            nucleotides downstream\
                                                            from the tested position\
                                                            as a training/testing\
                                                            feature',
                        required=False, default=20, dest="downstream")
    parser.add_argument('--tts-downstream', type=int, help='search for the TTS up to {--tts-downstream}\
                                                            from the gene 3\' end',
                        required=False, default=150, dest="tts_downstream")
    parser.add_argument('--tts-upstream', type=int, help='search for the TTS up to {--tts-upstream}\
                                                            from the gene 3\' end',
                        required=False, default=20, dest="tts_upstream")
    parser.add_argument('--tts-invasion', type=int, help='search for the TTS in the region with maximal overlap of {--tts-invasion} nt\
                                                            with the nearest downstream gene',
                        required=False, default=10, dest="tts_invasion")
    parser.add_argument('--min_probability', type=float, default=.15,
                        help='min probability of the classifier output to be considered a TTS',
                        required=False, dest="min_prob")
    parser.add_argument('--out-prefix', type=str, help='prefix of the output files',
                        required=False, default="termiRNAtor", dest="out_prefix")
    args = parser.parse_args() 

    return args