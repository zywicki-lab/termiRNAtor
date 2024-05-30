class ConfigFileHeaderError(Exception):
    """Raised when the header line in the config file is not properly
    formatted"""
    def __init__(self, line : str, message="Header line in the config file is not \
                                properly formatted. Should be \
                                'sample;path;condition;annotations;species;reference_genome;rna_seq_path;library_type'"):
        self.line = line
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.line} is not a proper header line." \
            "Should be 'sample;path;condition;annotations;species;reference_genome;rna_seq_path;library_type'"
            
class ConfigFileHeaderPredictError(Exception):
    """Raised when the header line in the config file is not properly
    formatted"""
    def __init__(self, line : str, message="Header line in the config file is not \
                                properly formatted. Should be \
                                'sample;annotations;species;reference_genome;rna_seq_path;library_type'"):
        self.line = line
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.line} is not a proper header line." \
            "Should be 'sample;annotations;species;reference_genome;rna_seq_path;library_type'"


class ConfigFileColumnError(Exception):
    """Raised when the line in the config cannot be properly divided '
    into eight columns"""

    def __init__(self, line : str, message="Line in the config file is not \
                                properly formatted. Should have eight columns \
                                separated by ';"):
        self.line = line
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.line} is not a properly formatted line." \
            "Should have seven columns separated by ';''"
            
class ConfigFilePredictionColumnError(Exception):
    """Raised when the line in the config cannot be properly divided '
    into eight columns"""

    def __init__(self, line : str, message="Line in the config file is not \
                                properly formatted. Should have six columns \
                                separated by ';"):
        self.line = line
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.line} is not a properly formatted line." \
            "Should have six columns separated by ';''"


class ConfigFileBAMNotFound(Exception):
    """Raised when the BAM file listed in the config coudln't be found"""

    def __init__(self, bam : str, message="BAM file listed in the config \
                                    coudln't be found"):
        self.bam = bam
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.bam} couldn't be found"


class ConfigFileGFFNotFound(Exception):
    """Raised when the GFF file listed in the config coudln't be found"""

    def __init__(self, gff : str, message="GFF file listed in the config \
                                    coudln't be found"):
        self.gff = gff
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.gff} couldn't be found"
    

class ConfigFileFASTANotFound(Exception):
    """Raised when the FASTA file listed in the config coudln't be found"""

    def __init__(self, fasta : str, message="FASTA file listed in the config \
                                    coudln't be found"):
        self.fasta = fasta
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.fasta} couldn't be found"


class ConfigFileInvalidBAMFile(Exception):
    """Raised when the BAM file is not properly formatted.
    Checked using samtools quickcheck by analysisng
    if the exit code is not equal to 0"""

    def __init__(self, bam : str, message="BAM file is not properly formatted"):
        self.bam = bam
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.bam} is not properly formatted"
    

class ConfigFileConditionNotNumeric(Exception):
    """Raised when in the condition column value other than numeric is provided"""

    def __init__(self, value : str, message="Condition value is invalid"):
        self.value = value
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.value} is not a proper condition value. Should be numeric"
    

class SpeciesNotFoundInConfigFile(Exception):
    """Raised when the species provided in the command line cannot be found in the ConfigFile"""

    def __init__(self, species : str, message="Species name not found in the ConfigFile"):
        self.species = species
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.species} not found in the ConfigFile" 
    
    
class ConfigFileWrongLibraryType(Exception):
    """Raised when the specified strand-specific library type is not R or F (single-end support only)"""

    def __init__(self, library : str, message="Library type must be one of ['R', 'F']"):
        self.library = library
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.library} - library type must be one of ['R', 'F']" 