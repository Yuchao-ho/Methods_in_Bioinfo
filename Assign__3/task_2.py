""" Author: Yuchao He 
    Aim of this file: To identify all 
    How to run the file:  python task_2.py
        """
import pandas as pd

input_data = pd.read_csv("ccm.csv")

# Identify the endpoints
endpoints = list(input_data[81:93]["To"].values)

def bfs(data:pd.DataFrame, ini_node:str) -> list:
    """
    Breadth-first search (BFS) to find all reachable metabolites from the initial metabolite.
    
    Parameters:
        data: A reduced reaction dataset.
        ini_node: The starting metabolite.

    Returns:
        list: A list of reachable metabolites.
    """
    p = []  # List to store visited nodes (metabolites)
    queue = [ini_node]  # Initialize the BFS queue with the initial node

    while len(queue) > 0:
        node = queue.pop(0)  
        queue.extend(data.loc[data["From"] == node]["To"].values)
        
        data2 = data.loc[(data["From"] != node) & (data["To"] != node)]
        data = data2
        # If the node hasn't been visited, add it to the visited list
        if node not in p:
            p.append(node)

    return p


def essential(data:pd.DataFrame, name:str, initial:str, goal:list) -> list: 
    """
    Find essential enzymes.

    Parameters:
        data: tested dataset.
        name (str): Name of the network.
        initial (str): The starting metabolite.
        goal (list): Endpoints.

    Outputs:
        Essential enzymes.
    """
    # Extract unique enzymes names in the network
    enzymes_0 = list(data["Enzyme"].values)
    enzymes = []
    for i in enzymes_0:
        if i not in enzymes:  
            enzymes.append(i)

    essen = []  # Store essential enzymes

    # Check each enzyme by using BSF
    for j in enzymes:
        # Remove all links related to enzyme j
        data2 = data.loc[data["Enzyme"] != j]
        p = bfs(data2, initial)
        for k in goal:
            if k not in p and j not in essen:
                essen.append(j)  

    print(f"The essential enzymes in {name} are:\n")
    for ele in essen:
        print(f"{ele}\n")


if __name__ == '__main__':
    essential(input_data, "central carbon metabolism", "Glucose [c]", endpoints)

