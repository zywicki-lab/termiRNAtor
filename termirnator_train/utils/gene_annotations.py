from typing import List, Dict, Tuple, Set
import gffutils
from . import configuration_file as cf


def create_gffutils_databases(samples : Dict[str, cf.ConfigSample]) -> Dict[str, gffutils.FeatureDB]:
    """ 

    Args:
        samples (Dict[str, ConfigSample]): samples dict returned by ConfigSample.parse_config

    Returns:
        Dict[str, gffutils.FeatureDB]: dict storing references to gffutils dfatabases
    """
    gff_files = [x.annotations for x in samples.values()]
    gff_files = set(gff_files)
    gff_databases = {}
    for gff in gff_files:
        # creates gffutils db in the same directory as a gff/gtf file
        db = gffutils.create_db(gff, dbfn=f'{gff}.db',
                                merge_strategy='merge', force=True,
                                keep_order=True, sort_attribute_values=True)
        gff_databases[f'{gff}.db'] = db
    return gff_databases