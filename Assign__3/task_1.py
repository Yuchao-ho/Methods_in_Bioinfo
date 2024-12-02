""" Author: Yuchao He 
    Aim of this file: To identify all non-essential enzymes in glycolysis.
    How to run the file:  python task_1.py
        """

import pandas as pd

## convert Fig_1 into .csv file,including start points, endpoints and related enzymes.
input_data = pd.read_csv("glycolysis.csv")

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


def non_essential(data, name, initial, goal):

    enzymes_0 = list(data["Enzyme"].values)
    enzymes = []
    for i in enzymes_0:
        if i not in enzymes:
            enzymes.append(i)

    nonessential = []
    for j in enzymes:
        data2 = data.loc[data["Enzyme"] != j]
        p = bfs(data2, initial)
        if goal in p:
            nonessential.append(j)

    print(f"The non-essential enzymes in {name} are:\n")
    for ele in nonessential:
        print(f"{ele}\n")


if __name__ == '__main__':
    non_essential(input_data, "glycolysis", "Glucose", "Pyruvate")