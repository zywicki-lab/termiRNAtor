
# termiRNAtor

<b>a tool for predicting intrinsic transcription termination signals based on bacterial RNA-seq data</b>

<hr>
<br>

termiRNAtor - predicts transcription termination sites in bacteria using machine learning approaches and RNA-seq data

<br>
<br>



## Introduction

Term-seq is a method designed for sequencing only 3’ ends of bacterial RNA molecules, which allows researchers to precisely identify transcription termination signals, RNA processing or decay products. 

Two variations of the Term-seq protocol have already been published. Detailed descriptions can be found in the following papers:

* Mondal S, Yakhnin AV, Sebastian A, Albert I, Babitzke P. NusA-dependent transcription termination prevents misregulation of global gene expression. Nat Microbiol. 2016 Jan 11;1:15007. doi: 10.1038/nmicrobiol.2015.7. PMID: 27571753; PMCID: PMC5358096.

* Dar D, Shamir M, Mellin JR, Koutero M, Stern-Ginossar N, Cossart P, Sorek R. Term-seq reveals abundant ribo-regulation of antibiotics resistance in bacteria. Science. 2016 Apr 8;352(6282):aad9822. doi: 10.1126/science.aad9822. PMID: 27120414; PMCID: PMC5756622.

The problem with the data originating from the Term-seq protocol is that only a minimal number of datasets and only from a highly limited set of species can be found in public databases. However, to truly understand the mechanisms of transcription termination across the bacterial kingdom, it is essential to analyze not only a few model species but a wide range of different organisms. This is why we came up with the idea to create a machine-learning model that can accurately predict transcription termination sites based solely on RNA-seq data. The creation of such a method was possible by analyzing Term-seq and RNA-seq samples originating from the same biological material for various bacterial species, incl. Bacillus subtilis, Escherichia coli, Listeria monocytogenes, Enterococcus faecalis, Zymomonas mobilis and Synechocystis sp. Term-seq data allows precise identification of the termination signals so that the appropriate ML model can be trained to predict them based only on the RNA-seq data and the genomic sequence.

termiRNAtor allows User to train their own RandomForest classifier based on paired Term-seq and RNA-seq data or use a pre-trained one (using the above-listed bacterial species).

The list of datasets employed in testing and training the classifier has been summarized in the file named `config.csv` (file names contain dataset IDs). This file contains information about both Term-seq and corresponding RNA-seq runs downloaded from the Sequence Read Archive. Data from Listeria monocytogenes have not been used for the purpose of classifier training but as independent dataset for model evaluation.

<br>
<br>
<hr>
<br>

## Training a new classifier - train_classifier.py

**NOTE: Use --help to learn more**

To train the new classifier, the User needs to specify for each sample:
- sorted BAM file describing mapped RNA-seq reads
- sorted BAM file describing mapped TERM-seq reads
- gene annotations in either GTF or GFF file formats (can be gzipped)
- reference genome in FASTA format (cannot be gzipped)
- RNA-seq library type (R - reverse-stranded, F - forward-stranded). Currently the tool supports only strand-specific single-end data
- species from which the data have originated (unique name, e.g. Staphylococcus_aureus)
- sample name or ID (unique)

Replicates belonging to the same sample should be assigned with a uniqe id (see `condition`). Only samples with at least two replicates can be used in the training process as the identification of the transcription termination sites is based on the signal consistency between the Term-seq replicates.

All above-described information should be encoded in the configuration file. An example can be found in the file named `config.csv`. Datasets present in the file were used to train the attached classifier.

Accuracy and other metrics of the pre-trained classifier can be found in:
- termiRNAtor_statistics.txt
- termiRNAtor_precission_recall.svg (precission and recall plot)
- termiRNAtor_precission_vs_recall.svg (precission vs recall plot)
- termiRNAtor_distances.svg (distribution of distances to the true TTS - see the algorithm)

<br>
<br>
<hr>

### Algorithm

1) Bam files are converted to bigwigs using `bamCoverage` (https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) from the `deeptools` package (--binSize 1 --normalizeUsing CPM). Separate coverage files are generated for each strand.
2) `TERMITe stable_rna_ends` and `annotate` (https://github.com/zywicki-lab/TERMITe-dev) are used to identify stable 3' RNA ends (only terminators predicted by both TransTermHP and RNAfold pipelines are taken into account)
3) To each gene, a region near its 3' end is assigned. The region starts by default 20 nt upstream from the 3' end (--tts-upstream) and ends 150nt downstream (--tts-downstream). However, overlaps with the following downstream gene located on the same strand are possible up to max 10 nt (--tts-invasion). In other case, the end will be equal to the start coordinate of the following CDS + --tts-invasion). Regions shorter than 50 nucleotides are excluded from the analysis.
4) The most abundant stable RNA end within the specified region is classified as a Transcription Termination Site (TTS). No TTS is assigned if no stable RNA end is reported for a given gene. Of course stable RNA ends must correspond to the TTS to be considered (predicted by both TransTermHP and RNAfold pipelines in TERMITe).
5) Genes with assigned TTS will be further used to train the RandomForestClassifier.
6) Both the sequence and the RNA-seq coverage (in 1nt bins) are extracted for each replicate from regions starting (--upstream) nt upstream (150 by default) and ending (--downstream) nt downstream (20 by default) from the identified TTS position. Each position within this range is used as a training example for machine learning algorithm. For each such position, we know from Term-seq whether it's classified as a TTS (see 4), therefore we can assign True/False label for each example. For each position, we also extract RNA-seq coverage and sequence as features to train the classifier. 
7) Genes having the average RNA-seq coverage of the first 20 nucleotides in the analyzed region (see 6) below median are excluded from the training and testing procedures as well as from the predictions.
8) Training set is created based on the data from bacterial species specified by the (--train_on), while the testing one is based on the (--test-on) argument. We suggest testing the model on species not included in the training procedure. By designing the experiment in this way, we know how the model behaves in species for which data were not included in the training process.
10) Positions marked as TTS in the training set are treated as positive examples (True label - see 6), while the rest as negative ones (False label - see 6). Both the positive and negative training sets are subsampled to have an equal size of (--training_set_size) training examples (if there are more than --training_set_size in a given set, subsampling with replacement is being used).
11) Testing and predictions are performed on the whole dataset from bacteria provided by the (--test-on) argument. The classifier is used to establish the probability of a position being a TTS for each gene region specified as in 3.). Only the position with the highest probability in the region and greater than (--min_probability) is treated as a TTS.

  * **TN (True Negative)** - number of genes without predicted TTS and without TTS established based on the Term-seq data
  * **FP (False Positive)** - number of genes for which predicted TTS is located > 5nt from the true TTS established based on Term-seq or with predicted TTS when no TTS has been established based on the Term-seq data
  * **TP (True Positive)** - number of genes for which predicted TTS is located <= 5nt from the true TTS established from Term-seq
  * **FN (False Negative)** - number of genes without predicted TTS but having the true TTS established based on Term-seq data

<br>
<br>
<hr>


### Example invocation
<br>

**NOTE: Currently project supports only single-end RNA-SEQ data**
<br>


`train_classifier.py -t config.csv -p 40 -u 150 -d 20 --train_on Enterococcus_faecalis,Escherichia_coli,Synechocystis,Zymomonas_mobilis,Bacillus_subtilis --test_on Listeria_monocytogenes`

Run
`train_classifier.py --help`
to learn more!


<br>
<br>
<hr>

### Results

* <b>termiRNAtor_classifier.pickle</b> - pickled RandomForestClassifier object
* <b>termiRNAtor_precission_vs_recall.svg</b> - precision vs recall plot
* <b>termiRNAtor_precission_recall.svg</b> - precision and recall vs probability
* <b>termiRNAtor_distances.svg</b> - histogram of distances between the predicted and true 3' RNA end
* <b>stdout</b> - classifier metrics and confusion matrix
* <b>termiRNAtor_statistics.txt</b> - same as information printed on stdout
<br>

`NOTE: you can replace "termiRNAtor" in the output file names with custom prefix using --out_prefix`

<hr>
<br>

## Use existing classifier - predict_TTS.py

**NOTE: Use --help to learn more**

if you want to use existing classifier to predict rho-independent termination sites based on RNA-seq data, use

`predict_TTS.py -c classifier.pickle -r config_prediction.csv -p 40 -u 150 -d 20`

## Funding

This work was supported by the Faculty of Biology Dean’s Grant (Adam Mickiewicz University in Poznan, Poland) `[GDWB-11/2020 to J.K.]`


