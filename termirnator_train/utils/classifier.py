from typing import List, Dict, Tuple, Set
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, precision_recall_curve
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle


def check_distance(row_with_max_prob, tts, return_distance = False):    
    gene = row_with_max_prob["gene"]
    sample = row_with_max_prob["rnaseq_sample"]
    try:
        for t in tts[f"{gene}x{sample}"]:
            if return_distance:
                return row_with_max_prob["genomic_position_pos"] - t
            if abs(row_with_max_prob["genomic_position_pos"] - t) <= 5:
                return True
    except KeyError:
        return False
    return False


def train_classifier(X_train : pd.DataFrame, X_labels : np.ndarray, n_estimators : int = 1000,
                    max_depth : int = None, criterion : str = "gini", min_samples_leaf : int = 1, class_weight : int = None,
                    threads : int = 1) -> RandomForestClassifier:
    """Train the RandomForestClassifier 
    
    Args:
        X_train (pd.DataFrame): training dataset
        X_labels (np.ndarray): training dataset labels
        n_estimators (int, optional): number of estimators to use. Defaults to 1000.
        max_depth (int, optional): max depth. Defaults to None.
        criterion (str, optional): criterion. Defaults to "gini".
        min_samples_leaf (int, optional): min samples leaf. Defaults to 1.
        class_weight (int, optional): class weight. Defaults to None.
        threads (int, optional): Number of threads. Defaults to None.
        

    Returns:
        RandomForestClassifier: classifier
    """
    
    X_labels = X_labels.astype('bool')
    classifier = RandomForestClassifier(n_estimators = n_estimators, random_state = 0, max_depth = max_depth, n_jobs = threads,
                                        criterion = criterion, min_samples_leaf = min_samples_leaf, class_weight = class_weight)
    classifier.fit(X_train, X_labels)

    return classifier


def draw_plots(classifier : RandomForestClassifier, Y_test : pd.DataFrame, Y_labels : pd.DataFrame,
                     Y_info : pd.DataFrame, Y_genes_to_keep : Dict[str, int], out_prefix : str):
    """Draw precission and recall vs probability

    Args:
        classifier (RandomForestClassifier): classifier
        Y_test (pd.DataFrame): testing dataset
        Y_labels (np.ndarray): testing labels
        Y_info (pd.DataFrame): testing dataset annotation
        Y_genes_to_keep (dict): genes to keep based on the coverage filtering applied in prepare_dataset
        out_prefix (str, optional): Prefix of the output files.
    """
    
    precision, recall = {}, {}
    predicted_proba = classifier.predict_proba(Y_test)
    for min_prob in np.arange(0.01, 1.0, 0.01):
        predictions = (predicted_proba [:,1] >= min_prob).astype('int')
        truth, labels, probabilities, distances = recalculate_labels(Y_test, Y_labels, Y_info, predictions, predicted_proba, Y_genes_to_keep, min_prob)
        a = confusion_matrix(truth, labels)
        tn, fp, fn, tp = a.ravel()
        try:
            p = (tp/(tp+fp))
            if np.isnan(p):
                precision[min_prob] = 1
            else:
                precision[min_prob] =  p
        except ZeroDivisionError:
            precision[min_prob] = 1
        try:
            r = tp/(tp+fn)
            recall[min_prob] = r
        except ZeroDivisionError:
            recall[min_prob] = np.nan

    plt.figure(figsize=(5, 2.7), layout='constrained')
    plt.plot(precision.keys(), precision.values(), label='precision')
    plt.plot(recall.keys(), recall.values(), label='recall')
    plt.xlabel('probability')
    plt.ylabel('value')
    plt.title("Precission and recall")
    plt.legend()
    plt.savefig(f"{out_prefix}_precission_recall.svg")
    
    plt.figure(layout='constrained')
    plt.plot(list(precision.values()), list(recall.values()), label='RF classifier')
    plt.plot([0,1], [1,0], "r--", label='random')
    plt.xlabel('Precision')
    plt.ylabel('Recall')
    plt.xlim([min(precision.values()),max(precision.values())])
    plt.title("Precission vs recall")
    plt.legend()
    plt.savefig(f"{out_prefix}_precission_vs_recall.svg")
 
        
def draw_distances_histogram(distances : List[int], out_prefix : str) -> None:
    """Draws histogram of distances between the predicted and true 3' RNA end

    Args:
        distances (List[int]): list of distances
        out_prefix (str, optional): Prefix of the output files.
    """
    fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
    n, bins, patches = ax.hist(distances, 300, density=True, facecolor='C0', alpha=0.75)
        
    ax.set_xlabel('distance [nt]')
    ax.set_ylabel('Probability')
    ax.set_title("Distance to the true 3' end")
    ax.set_xlim(-20, 20)
    ax.grid(True)
    plt.savefig(f"{out_prefix}_distances.svg")


def recalculate_labels(Y_test : pd.DataFrame, Y_labels : np.ndarray, test_info : pd.DataFrame, predictions : np.ndarray, proba : np.ndarray,
                       genes_to_keep_test : dict, min_probability : float) -> tuple:
    """Recalculate the labels so as a true predictions are taken those above probability of min_probability

    Args:
        Y_test (pd.DataFrame): testing dataset
        Y_labels (np.ndarray): testing labels
        test_info (pd.DataFrame): testing dataset annotation
        predictions (np.ndarray): model predictions
        proba (np.ndarray): probabilities of beeing TTS
        genes_to_keep_test (dict): genes to keep based on the coverage filtering applied in prepare_dataset
        min_probability (float): probability threshold for decision making

    Returns:
        tuple: testing dataset, recalculated labels, additional_information about the classifier, max probabilities and distances from the predicted to true TTS
    """
    
    Z_test, Z_labels, tts = [], [], {}
    
    test_dataset = pd.DataFrame(Y_test, columns=[f"col{x}" for x in range(Y_test.shape[1])] )
    test_dataset["prediction"] = predictions
    test_dataset["proba"] = proba[:,1]
    test_dataset["truth"] = Y_labels
    test_dataset = pd.merge(test_info, test_dataset, left_index=True, right_index=True)
    del test_dataset["cov_0"]
    
    genes = list(set(test_dataset["gene"]))  # list of all genes
    
    ### TTS dictionary
    positive_test_samples = test_dataset[test_dataset["truth"] == True]
    for index, row in positive_test_samples.iterrows():    
        if f"{row['gene']}x{row['rnaseq_sample']}" not in tts:
            tts[f"{row['gene']}x{row['rnaseq_sample']}"] = [row['genomic_position_pos']]   
        else:
            tts[f"{row['gene']}x{row['rnaseq_sample']}"].append([row['genomic_position_pos']])
    
    os.system("mkdir -p predictions")
    os.system("rm predictions/*")
    probabilities = []
    distances = []
    for gene in genes:
        gene_df = test_dataset.loc[test_dataset['gene'] == gene]
        samples = list(set(gene_df["rnaseq_sample"]))  # list of all samples
        for sample in samples:
            if f"{gene}{sample}" not in genes_to_keep_test:
                continue
            #with open(f"predictions/{sample.split('/')[-1].replace('.bam', '')}.bed", "a") as f_out:  # TODO add tab below if to be reinstated
            sample_df = gene_df.loc[gene_df['rnaseq_sample'] == sample]
            row_with_max_prob = sample_df.loc[sample_df['proba'].idxmax()]
            probabilities.append(row_with_max_prob["proba"])
            if row_with_max_prob["proba"] >= min_probability:
                distances.append(check_distance(row_with_max_prob, tts, True))
                if check_distance(row_with_max_prob, tts):
                    #print(f"{row_with_max_prob['genomic_position_chr']}\t{row_with_max_prob['genomic_position_pos']}\t{row_with_max_prob['genomic_position_pos']}\t{gene}_TP\t{row_with_max_prob['proba']}\t{row_with_max_prob['genomic_position_strand']}", file=f_out)
                    Z_test.append(True)
                    Z_labels.append(True)
                else:
                    #print(f"{row_with_max_prob['genomic_position_chr']}\t{row_with_max_prob['genomic_position_pos']}\t{row_with_max_prob['genomic_position_pos']}\t{gene}_FP\t{row_with_max_prob['proba']}\t{row_with_max_prob['genomic_position_strand']}", file=f_out)
                    Z_test.append(False)
                    Z_labels.append(True)
            elif row_with_max_prob["proba"] < min_probability and True not in list(sample_df.truth):
                #print(f"{row_with_max_prob['genomic_position_chr']}\t{row_with_max_prob['genomic_position_pos']}\t{row_with_max_prob['genomic_position_pos']}\t{gene}_TN\t{row_with_max_prob['proba']}\t{row_with_max_prob['genomic_position_strand']}", file=f_out)
                Z_test.append(False)
                Z_labels.append(False)
            else:
                #print(f"{row_with_max_prob['genomic_position_chr']}\t{row_with_max_prob['genomic_position_pos']}\t{row_with_max_prob['genomic_position_pos']}\t{gene}_FN\t{row_with_max_prob['proba']}\t{row_with_max_prob['genomic_position_strand']}", file=f_out)
                Z_test.append(True)
                Z_labels.append(False)   
    return Z_test, Z_labels, probabilities, distances


def run_classifier_training(X_train : pd.DataFrame, X_labels : pd.DataFrame, Y_test : pd.DataFrame, Y_labels : pd.DataFrame,
                     Y_info : pd.DataFrame, Y_genes_to_keep : Dict[str, int], out_prefix : str, min_probability : int = 0.15,  threads : int = 1,
                     ):
    """Trains the RandomForestClassifier and evaluates the model efficiency

    Args:
        X_train (pd.DataFrame): training dataset
        X_labels (np.ndarray): training labels
        Y_test (pd.DataFrame): testing dataset
        Y_labels (np.ndarray): testing labels
        Y_info (pd.DataFrame): testing dataset annotation
        Y_genes_to_keep (Dict[str, int]): dictionary of the genes to be kept based on their expression leveles (output of the filter_by_avg_exp)
        min_probability (int, optional): min probability of the classifier output to be considered TTS. Defaults to 0.15.
        threads (int, optional): number of threads to use. Defaults to 1.
        out_prefix (str, optional): Prefix of the output files.
        

    Returns:
        RandomForestClassifier: classifier object
    """
    
    classifier = train_classifier(X_train, X_labels, threads = threads)
              
    predicted_proba = classifier.predict_proba(Y_test)  # predict probabilities of beeing TTS for each testing example

    draw_plots(classifier, Y_test, Y_labels, Y_info, Y_genes_to_keep, out_prefix)
    predictions = (predicted_proba [:,1] >= min_probability).astype('int')
    truth, labels, probabilities, distances = recalculate_labels(Y_test, Y_labels, Y_info, predictions, predicted_proba, Y_genes_to_keep, min_probability)
    draw_distances_histogram(distances, out_prefix)
    print(confusion_matrix(truth, labels))
    tn, fp, fn, tp = confusion_matrix(truth, labels).ravel()
    print(f"TP: {tp}; TN: {tn}; FN: {fn}; FP: {fp}")
    print(classification_report(truth, labels))
    
    # save output
    with open(f'{out_prefix}_classifier.pickle', 'wb') as handle:
        pickle.dump(classifier, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    with open(f"{out_prefix}_statistics.txt", "w") as f:
        print(confusion_matrix(truth, labels), file=f)
        print(f"TP: {tp}; TN: {tn}; FN: {fn}; FP: {fp}", file=f)
        print(classification_report(truth, labels), file=f)
    
    return classifier


def predict_using_classiffer(classifier : RandomForestClassifier, info : pd.DataFrame, dataset : pd.DataFrame,
                             genes_to_keep : Dict[str, float], out_prefix : str) -> None:
    """Uses classifier to predict TTS based on the RNA-seq data

    Args:
        classifier (RandomForestClassifier): random forest classifier
        info (pd.DataFrame): dataframe with gene info corresponding to the dataset rows
        dataset (pd.DataFrame): dataset
        genes_to_keep (Dict[str, float]): genes with sufficient expression levels to be kept 
    """
    predicted_proba = classifier.predict_proba(dataset)
    
    dataset_merged = pd.merge(info, dataset, left_index=True, right_index=True)
    
    os.system(f"mkdir -p {out_prefix}_predictions")
        
    del dataset_merged["cov_0_x"]
    dataset_merged = dataset_merged.rename(columns={'cov_0_y': 'cov_0'})
    dataset_merged["proba"] = predicted_proba[:,1]
    genes = list(set(dataset_merged["gene"]))
    samples = list(set(dataset_merged["rnaseq_sample"]))
    d_samples = {}

    for sample in samples:
        s = sample.split("/")[-1].replace(".bam", "") + ".bed"
        d_samples[sample] = open(f"{out_prefix}_predictions/{s}", "w")

    for gene in genes:
        gene_df = dataset_merged.loc[dataset_merged['gene'] == gene]
        samples = list(set(gene_df["rnaseq_sample"]))
        for sample in samples:
            if f"{gene}{sample}" not in genes_to_keep:
                continue
            sample_df = gene_df.loc[gene_df['rnaseq_sample'] == sample]
            row_with_max_prob = sample_df.loc[sample_df['proba'].idxmax()]
            if row_with_max_prob["proba"] >= 0.1:
                chromosome = row_with_max_prob["genomic_position_chr"]
                position = row_with_max_prob["genomic_position_pos"]
                strand = row_with_max_prob["genomic_position_strand"]
                name = row_with_max_prob["gene"]
                score = row_with_max_prob["proba"]
                print(f"{chromosome}\t{position}\t{position}\t{name}\t{score}\t{strand}", file=d_samples[sample])