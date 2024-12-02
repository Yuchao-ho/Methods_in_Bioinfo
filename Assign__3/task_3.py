""" Author: Yuchao He 
    Aim of this file: To identify all 
    How to run the file:  python task_3.py
        """

import pandas as pd

def get_data(data_path, rna_path):
    """ 
    Load and preprocess the metabolic network and RNA count data. 
    """

    input_data = pd.read_csv(data_path)

    rnacounts = pd.read_csv(rna_path)

    healthy_data = pd.merge(input_data, rnacounts[["Enzyme", "Healthy RNA count"]], on="Enzyme", how="left")
    healthy_data = healthy_data.loc[healthy_data["Healthy RNA count"] != 0]

    cancer_data = pd.merge(input_data, rnacounts[["Enzyme", "Cancer RNA count"]], on="Enzyme", how="left")
    cancer_data = cancer_data.loc[cancer_data["Cancer RNA count"] != 0]

    endpoints = list(input_data[81:93]["To"].values)

    return healthy_data, cancer_data, endpoints


def bfs(data, initial):

    p = []
    queue = [initial]
    while len(queue) > 0:
        node = queue.pop(0)
        queue.extend(data.loc[data["From"] == node]["To"].values)
        data2 = data.loc[(data["From"] != node) & (data["To"] != node)]
        data = data2
        if node not in p:
            p.append(node)

    return p


def target_enzy(data_path, rna_path, initial):
    """ 
    Identify enzymes for targeting cancer cells. 
    """

    healthy_data, cancer_data, endpoints = get_data(data_path, rna_path)

    healthy_enzymes_0 = list(healthy_data["Enzyme"].values)
    healthy_enzymes = []
    for i in healthy_enzymes_0:
        if i not in healthy_enzymes:
            healthy_enzymes.append(i)

    cancer_enzymes_0 = list(cancer_data["Enzyme"].values)
    cancer_enzymes = []
    for i in cancer_enzymes_0:
        if i not in cancer_enzymes:
            cancer_enzymes.append(i)

    # Identify essential enzymes for cancer cells
    cancer_essential = []
    for j in cancer_enzymes:
        cancerdata2 = cancer_data.loc[cancer_data["Enzyme"] != j]
        p = bfs(cancerdata2, initial)
        for k in endpoints:
            if k not in p and j not in cancer_essential:
                cancer_essential.append(j)

    # Identify non-essential enzymes for healthy cells
    healthy_nonessential = []
    for j in healthy_enzymes:
        healthydata2 = healthy_data.loc[healthy_data["Enzyme"] != j]
        p = bfs(healthydata2, initial)
        if set(endpoints).issubset(p) and j not in healthy_nonessential:
            healthy_nonessential.append(j)

    # Target enzymes: essential in cancer cells & non-essential in healthy cells
    solution = [
        enzyme
        for enzyme in cancer_essential
        if enzyme in healthy_nonessential or enzyme not in healthy_enzymes
    ]

    print(f"The enzymes suitable to kill cancer-cells are:\n")
    for ele in solution:
        print(f"{ele}\n")
    


if __name__ == '__main__':
    data_path = "ccm.csv"
    rna_path = "RNAcounts.csv"
    target_enzy(data_path, rna_path, "Glucose [c]")