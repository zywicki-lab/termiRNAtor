from . import exceptions
from os.path import isfile, exists
from subprocess import Popen, PIPE
from io import TextIOWrapper
import os
from typing import Union, Dict, Tuple, List
from collections import Counter
import itertools
from random import choices
from . import miscellaneous
import string
from copy import deepcopy


def parse_termirnator_path(input : Union[str, list]) -> Union[str, list]:
    """Replace $TERMIRNATOR in the config file with a proper path to the termirnator project
    
    Args:
        input (Union): single line or list of lines

    Returns:
        Union: str or list with $TERMIRNATOR replaced with a proper path to the termirnator dir
    """
    path = os.path.abspath(__file__)
    utils_path = os.path.abspath(os.path.join(path, os.pardir))
    termirnator_train_path = os.path.abspath(os.path.join(utils_path, os.pardir))
    termirnator_path = os.path.abspath(os.path.join(termirnator_train_path, os.pardir))
    
    if isinstance(input, str):
        # if the input argument represents a single line from the config file
        if "$TERMIRNATOR" in input:
            input = input.replace("$TERMIRNATOR", termirnator_path)
    elif isinstance(input, list):
        for line_idx in range(0, len(input)):
            if "$TERMIRNATOR" in input[line_idx]:
                input[line_idx] = input[line_idx].replace("$TERMIRNATOR", termirnator_path)
                
    return input
    

class ConfigSample():
    def __init__(self, sample_id, bam_path, condition, annotations, species,
                 ref_genome, rna_seq_path, library):
        self.sample_id = sample_id  # ID
        self.bam_path = bam_path  # path to BAM file
        self.condition = condition  # replicate group
        self.annotations = annotations  # path to the GTF/GFF3
        self.species = species  # species
        self.ref_genome = ref_genome  # reference genome in the fasta file format
        self.rna_seq_path = rna_seq_path  # path to the corresponding RNA-seq BAM file
        self.library = library # RNA-seq library type (either F or R)
        
    def __str__(self):
        return f"{self.sample_id}"
    
    def __repr__(self):
        return f"{self.sample_id}"
    
    def print_line(self):
        return f"{self.sample_id};{self.bam_path};{self.condition};{self.annotations};{self.species};{self.ref_genome};{self.ref_genome};{self.rna_seq_path};{self.library}"

    
    @staticmethod
    def bam_quickcheck(bam : str):
        """Runs samtools quickcheck to validate BAM

        Args:
            bam (str): path to bam file

        Raises:
            exceptions.ConfigFileInvalidBAMFile: [description]
        """
        p = Popen(["samtools", "quickcheck", f"{bam}"],
                      stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            raise exceptions.ConfigFileInvalidBAMFile(bam)
    
    @staticmethod
    def validate_config_file(config_file : TextIOWrapper):
        # check the header
        header = config_file.readline()
        if header.strip() != "sample;path;condition;annotations;species;reference_genome;rna_seq_path;library_type":
            raise exceptions.ConfigFileHeaderError(header)
        lines = config_file.readlines()
        lines = parse_termirnator_path(lines)
   
        for line in lines:
            if line.startswith("#"):
                continue
            sl = line.strip().split(";")
            if len(sl) != 8: # check the number of columns
                raise exceptions.ConfigFileColumnError(line)
            if not isfile(sl[1]):
                raise exceptions.ConfigFileBAMNotFound(sl[1])
            if not isfile(sl[3]):
                raise exceptions.ConfigFileGFFNotFound(sl[3])
            if not isfile(sl[5]):
                raise exceptions.ConfigFileFASTANotFound(sl[5])
            if not isfile(sl[6]):
                raise exceptions.ConfigFileBAMNotFound(sl[6])
            ConfigSample.bam_quickcheck(sl[1])
            ConfigSample.bam_quickcheck(sl[6])
            if not sl[2].isnumeric():
                raise exceptions.ConfigFileConditionNotNumeric(sl[2])
            if sl[7] not in ["R", "F"]:
                raise exceptions.ConfigFileWrongLibraryType(sl[7])
    
    @staticmethod
    def is_header_or_comment(line : str) -> bool:
        """Returns True if the line from the config file is a header line or comment

        Args:
            line (str): line

        Returns:
            bool: True if header or comment, False otherwise
        """
        if line.startswith("sample;") or line.startswith("#"):
            return True
        else:
            return False
        
    @staticmethod
    def parse_line(line : str) -> Tuple[str, object]:
        """Parse a single line from the config file and returns ConfigSample object

        Args:
            line (str): data line from the config file

        Returns:
            Tuple[str, ConfigSample]: tuple with sample id and ConfigSample object;
                                      None if header or comment line or failed to parse
        """
        if ConfigSample.is_header_or_comment(line):
            return None
        try:
            sl = line.strip().split(";")
            sample_id = sl[0]  # ID
            bam_path = sl[1]  # path to BAM file
            condition = sl[2]  # replicate group
            annotations = sl[3]  # path to the GTF/GFF3
            species = sl[4]  # species
            ref_genome = sl[5]  # reference genome
            rna_seq_path = sl[6]  # path to the corresponding RNA-seq BAM
            library = sl[7]  # library type
            return_tuple = ( sample_id, ConfigSample(sample_id, bam_path, condition,
                                            annotations, species, ref_genome,
                                            rna_seq_path, library) )
        except IndexError:
            miscellaneous.eprint("Invalid format of the config file")
            return None
        return return_tuple
        
            
    @staticmethod
    def parse_config(config_file : TextIOWrapper) -> Dict[str, object]:
        """
        Parse configuration file describing Term-seq sample.
        An example of a properly written config file can be found
        in sample_data/config.csv.
        TERMITe should be run with exactly 2 replicates per sample.
        If more replicates are present, each pairwise combination
        of them will be included in the modified config data as a separate condition

        Args:
            config_file (TextIOWrapper): configuration file

        Returns:
            Dict[str, object]: list of the modified config lines
        """
        
        ConfigSample.validate_config_file(config_file)
        config_file.seek(0)  # get back to the first byte of the file
        
        config = {}
        
        for line in config_file:
            if not ConfigSample.is_header_or_comment(line):   # not the header line nor the comment line
                sample = ConfigSample.parse_line(line)
                config[sample[0]] = sample[1]
                
        conditions = [sample.condition for sample_id, sample in config.items()]
        counter = Counter(conditions)
        conditions_to_split = [key for key, value in counter.items() if value > 2]
        conditions_not_affected = {sample_id : sample for sample_id, sample in config.items() if sample.condition not in conditions_to_split}
        conditions_splitted = {}
        for condition in conditions_to_split:
            replicates = [sample for sample in config.values() if sample.condition == condition]
            indexes = list(range(len(replicates)))
            combinations = itertools.combinations(indexes, 2)
            for combination in combinations:
                random_replicate_suffix = "".join(choices(string.ascii_lowercase, k=5))
                while f"{condition}_{random_replicate_suffix}" in conditions:
                    random_replicate_suffix = "".join(choices(string.ascii_lowercase, k=5))
                mod_repl_1 = deepcopy( replicates[combination[0]] )
                mod_repl_1.sample_id = f"{mod_repl_1.sample_id}_{random_replicate_suffix}"
                mod_repl_1.condition = f"{condition}_{random_replicate_suffix}"
                mod_repl_2 = deepcopy( replicates[combination[1]] )
                mod_repl_2.condition = f"{condition}_{random_replicate_suffix}"
                mod_repl_2.sample_id = f"{mod_repl_2.sample_id}_{random_replicate_suffix}"
                conditions_splitted[mod_repl_1.sample_id] = mod_repl_1
                conditions_splitted[mod_repl_2.sample_id] = mod_repl_2
        
        conditions_splitted.update(conditions_not_affected)
            
        modified_config = conditions_splitted
        
        return modified_config
    
    
class ConfigSamplePredict(ConfigSample):
    def print_line(self):
        return f"{self.sample_id};{self.condition};{self.annotations};{self.species};{self.ref_genome};{self.ref_genome};{self.rna_seq_path};{self.library}"

    @staticmethod
    def validate_config_file(config_file : TextIOWrapper):
        # check the header
        header = config_file.readline()
        if header.strip() != "sample;annotations;species;reference_genome;rna_seq_path;library_type":
            raise exceptions.ConfigFileHeaderPredictError(header)
        lines = config_file.readlines()
        lines = parse_termirnator_path(lines)
   
        for line in lines:
            if line.startswith("#"):
                continue
            sl = line.strip().split(";")
            if len(sl) != 6: # check the number of columns
                raise exceptions.ConfigFilePredictionColumnError(line)
            if not isfile(sl[1]):
                raise exceptions.ConfigFileGFFNotFound(sl[2])
            if not isfile(sl[3]):
                raise exceptions.ConfigFileFASTANotFound(sl[3])
            if not isfile(sl[4]):
                raise exceptions.ConfigFileBAMNotFound(sl[4])
            ConfigSamplePredict.bam_quickcheck(sl[4])
            if sl[5] not in ["R", "F"]:
                raise exceptions.ConfigFileWrongLibraryType(sl[5])
            
    @staticmethod
    def parse_line(line : str) -> Tuple[str, object]:
        """Parse a single line from the config file and returns ConfigSamplePredict object

        Args:
            line (str): data line from the config file

        Returns:
            Tuple[str, ConfigSamplePredict]: tuple with sample id and ConfigSamplePredict object;
                                      None if header or comment line or failed to parse
        """
        if ConfigSamplePredict.is_header_or_comment(line):
            return None
        try:
            sl = line.strip().split(";")
            sample_id = sl[0]  # ID
            condition = None
            annotations = sl[1]  # path to the GTF/GFF3
            species = sl[2]  # species
            ref_genome = sl[3]  # reference genome
            rna_seq_path = sl[4]  # path to the corresponding RNA-seq BAM
            library = sl[5]  # library type
            return_tuple = ( sample_id, ConfigSamplePredict(sample_id, None, condition,
                                            annotations, species, ref_genome,
                                            rna_seq_path, library) )
        except IndexError:
            miscellaneous.eprint("Invalid format of the config file")
            return None
        return return_tuple
        
            
    @staticmethod
    def parse_config(config_file : TextIOWrapper) -> Dict[str, object]:
        """
        Parse configuration file describing RNA-seq sample.
        An example of a properly written config file can be found
        in config_prediction.csv.

        Args:
            config_file (TextIOWrapper): configuration file

        Returns:
            Dict[str, object]: list of the modified config lines
        """
        
        ConfigSamplePredict.validate_config_file(config_file)
        config_file.seek(0)  # get back to the first byte of the file
        
        config = {}
        
        for line in config_file:
            if not ConfigSamplePredict.is_header_or_comment(line):   # not the header line nor the comment line
                sample = ConfigSamplePredict.parse_line(line)
                config[sample[0]] = sample[1]
    
        return config

