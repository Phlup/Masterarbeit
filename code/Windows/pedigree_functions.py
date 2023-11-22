import tskit
import msprime
import stdpopsim
import networkx as nx
import matplotlib.pyplot as plt
import io
from IPython.display import display as SVG
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict

#generate pedigree df with n "founder" individuals
def pedigree_init(n: int = 10) -> pd.DataFrame:
    """
    Generate a DataFrame representing a pedigree with 'n' founder individuals.

    Parameters:
        n (int): Number of founder individuals.

    Returns:
        pd.DataFrame: DataFrame representing the pedigree.
    """
    pedigree = pd.DataFrame(
    {"id": [i for i in range(0, n)],
    "parent0": ["."]*n,
    "parent1": ["."]*n,
    "time": [0]*n})
    return pedigree

#selfing (adding selfed individual to pedigree)
def add_selfing(df: pd.DataFrame, size_diff: int = 0) -> pd.DataFrame:
    """
    Add selfed individuals (next generation) to the pedigree.

    Parameters:
        df (pd.DataFrame): DataFrame representing the pedigree.
        size_diff (int): Size difference of current and next generation.

    Returns:
        pd.DataFrame: Updated DataFrame with selfed individuals.
    """
    df['time'] += 1
    parents = df.loc[df['time'] == df['time'].min(), "id"].to_numpy()
    max_id = df.loc[df['id'].max(), "id"]
    ids = [i for i in range(max_id + 1, max_id + 1 + len(parents) + size_diff)]
    df = pd.concat(objs=(df, pd.DataFrame({"id" : ids, "parent0": np.resize(parents, len(parents) + size_diff),
                                           "parent1": np.resize(parents, len(parents) + size_diff),
                                           "time" : [0]*(len(parents) + size_diff)}))).reset_index(drop = True)
    return(df)

#random mating
#implement chance of selfing
def add_random_mating(df: pd.DataFrame, size_diff: int = 0) -> pd.DataFrame:
    """
    Add individuals from random mating (next generation) to the pedigree.

    Parameters:
        df (pd.DataFrame): DataFrame representing the pedigree.
        size_diff (int): Size difference of current and next generation.

    Returns:
        pd.DataFrame: Updated DataFrame with randomly mated individuals.
    """
    df['time'] += 1
    parents = df.loc[df['time'] == df['time'].min(), "id"].to_numpy()
    max_id = df.loc[df['id'].max(), "id"]
    ids = [i for i in range(max_id + 1, max_id + 1 + len(parents) + size_diff)]
    df = pd.concat(objs=(df, pd.DataFrame({"id" : ids, "parent0": np.random.choice(parents, len(parents) + size_diff),
                                           "parent1": np.random.choice(parents, len(parents) + size_diff),
                                           "time" : [0]*(len(parents) + size_diff)}))).reset_index(drop = True)
    
    return(df)
    
#selective mating
def add_selective_mating(df: pd.DataFrame, parents: Tuple[int, int], offspring: int = 1) -> pd.DataFrame:
    """
    Add selectively mated individuals to the pedigree.

    Parameters:
        df (pd.DataFrame): DataFrame representing the pedigree.
        parents (Tuple[int, int]): IDs of the parents.
        offspring (int): Number of offspring.

    Returns:
        pd.DataFrame: Updated DataFrame with selectively mated individuals.
    """
    df['time'] += 1
    max_id = df.loc[df['id'].max(), "id"]
    ids = [i for i in range(max_id + 1, max_id + 1 + offspring)]
    df = pd.concat(objs=(df, pd.DataFrame({"id": ids, "parent0": np.resize(parents[0], offspring),
                                           "parent1": np.resize(parents[1], offspring),
                                           "time": [0]*offspring}))).reset_index(drop = True)     
    return(df)

def cross_selfing_ped(offspring: int, selfing_genos: int = 5) -> pd.DataFrame:
    """
    Set up pedigree for crossing two founders with arbitrary offspring and selfing generations.

    Parameters:
        offspring (int): Number of offspring.
        selfing_genos (int): Number of selfing generations.

    Returns:
        pd.DataFrame: DataFrame representing the pedigree.
    """
    cs_ped_df = pedigree_init(n = 2) 
    cs_ped_df = add_selective_mating(cs_ped_df, list([0,1]), offspring = offspring)
    for i in range(0, selfing_genos):
        cs_ped_df = add_selfing(cs_ped_df)
    return(cs_ped_df)

#generate msprime trees from df
def df_to_ts(df: pd.DataFrame, seq_len: int = 100) -> msprime.TreeSequence:
    """
    Generate a tree sequence from a DataFrame representing a pedigree.

    Parameters:
        df (pd.DataFrame): DataFrame representing the pedigree.
        seq_len (int): Length of the sequence/chromosome.

    Returns:
        msprime.TreeSequence: Generated tree sequence.
    """
    ped_string = df.to_string(index = False)
    ts_ped = msprime.parse_pedigree(io.StringIO("#" + ped_string), sequence_length = seq_len)
    return ts_ped

