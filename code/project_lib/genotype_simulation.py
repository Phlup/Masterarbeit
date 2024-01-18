import tskit
import msprime
import stdpopsim
import networkx as nx
import matplotlib.pyplot as plt
import io
from IPython.display import display as SVG
from decimal import *
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict
import scipy.stats as stats
import seaborn as sns

##pedigree functions
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
def add_random_mating(df: pd.DataFrame, size_diff: int = 0, selfing_rate: int = 0) -> pd.DataFrame:
    """
    Add individuals from random mating (next generation) to the pedigree.

    Parameters:
        df (pd.DataFrame): DataFrame representing the pedigree.
        size_diff (int): Size difference of current and next generation.
        selfing_rate (int): Specifies rate of guaranteed selfing.

    Returns:
        pd.DataFrame: Updated DataFrame with randomly mated individuals.
    """
    df['time'] += 1
    parents = df.loc[df['time'] == df['time'].min(), "id"].to_numpy()
    max_id = df.loc[df['id'].max(), "id"]
    ids = [i for i in range(max_id + 1, max_id + 1 + len(parents) + size_diff)]
    selfs = int(Decimal((len(parents) + size_diff)*selfing_rate).quantize(Decimal('1'), rounding=ROUND_HALF_UP))
    non_selfs = int(Decimal((len(parents) + size_diff)*(1-selfing_rate)).quantize(Decimal('1'), rounding=ROUND_HALF_DOWN))
    parent0 = np.concatenate([np.resize(parents, selfs), np.random.choice(parents, non_selfs)])
    parent1 = np.concatenate([np.resize(parents, selfs), np.random.choice(parents, non_selfs)])
    df = pd.concat(objs=(df, pd.DataFrame({"id" : ids, "parent0": parent0,
                                           "parent1": parent1,
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
def df_to_ts(df: pd.DataFrame, seq_len: int = 100) -> tskit.TreeSequence:
    """
    Generate a tree sequence from a DataFrame representing a pedigree.

    Parameters:
        df (pd.DataFrame): DataFrame representing the pedigree.
        seq_len (int): Length of the sequence/chromosome.

    Returns:
        tskit.TreeSequence: Generated tree sequence.
    """
    ped_string = df.to_string(index = False)
    ts_ped = msprime.parse_pedigree(io.StringIO("#" + ped_string), sequence_length = seq_len)
    return ts_ped

##genotype simulation functions
#propagate geno functions
def get_set(ts: tskit.TreeSequence) -> List[int]:
    """
    Get a list of all unique edges in the tree sequence (each edge has at least one child or parent).

    Parameters:
        ts (tskit.TreeSequence): A tree sequence.

    Returns:
        List[int]: List of unique edges.
    """
    return(list(set(ts.edges_parent).union(set(ts.edges_child))))

def get_founders(ts: tskit.TreeSequence) -> List[int]:
    """
    Get a list of founders in the tree sequence (all nodes except founders have parents).

    Parameters:
        ts (tskit.TreeSequence): A tree sequence.

    Returns:
        List[int]: List of founder nodes.
    """
    return(list(set(ts.edges_parent) - set(ts.edges_child)))

def get_offspring(ts: tskit.TreeSequence) -> List[int]:
    """
    Get a list of offspring (last generation) in the tree sequence (all nodes except offspring have children).

    Parameters:
        ts (tskit.TreeSequence): A tree sequence.

    Returns:
        List[int]: List of offspring nodes.
    """
    return(list(set(ts.edges_child) - set(ts.edges_parent)))

def get_edges(parents: List[int], ts: tskit.TreeSequence) -> pd.DataFrame:
    """
    Get a DataFrame containing all edges of next generation in ARG.

    Parameters:
        parents (List[int]): List of parent nodes.
        ts (tskit.TreeSequence): A tree sequence.

    Returns:
        pd.DataFrame: DataFrame containing edges.
    """
    pc_df = pd.DataFrame({"left": ts.edges_left, "right": ts.edges_right, "parent" : ts.edges_parent, "child" : ts.edges_child})
    children = pc_df.loc[pc_df["parent"].isin(parents),]
    return(children)

#function to choose founders and split diploid genome into ts nodes -> output genoypes df to use in propagate_geno function
#not generalizable, specific to dataset from bukler et al 2007
def get_founder_nodes(genotypes: pd.DataFrame, founder_list: List[str]) -> pd.DataFrame:
    """
    Get founder nodes from genotypes DataFrame.

    Parameters:
        genotypes (pd.DataFrame): DataFrame containing genotypes.
        founder_list (List[str]): List of two founder names.
founders
    Returns:
        pd.DataFrame: DataFrame containing founder nodes.
    """
    founder_geno = genotypes.loc[genotypes['RIL'].isin(founder_list)].drop("RIL", axis = 1).reset_index(drop = True)
    founder_nodes = pd.concat([founder_geno.loc[0].str[0],
                              founder_geno.loc[0].str[2],
                              founder_geno.loc[1].str[0],
                              founder_geno.loc[1].str[2]], axis = 1).T.reset_index(drop = True)
    return founder_nodes

#propagate genotypes of founders along ancestral recombination graph using genetic map
def propagate_geno(arg: tskit.TreeSequence, founder_nodes: pd.DataFrame, genmap: pd.DataFrame) -> pd.DataFrame:
    """
    Propagate genotypes through the tree sequence/ARG.

    Parameters:
        arg (tskit.TreeSequence): TreeSequence object obtained from msprime.sim_ancestry.
        founder_nodes (pd.DataFrame): DataFrame containing founder nodes.
        genmap (pd.DataFrame): DataFrame containing genetic map.

    Returns:
        pd.DataFrame: DataFrame containing propagated genotypes.
    """
    sites = genmap["Marker"]
    geno_sim = pd.DataFrame(columns = list(["node"]) + list(sites))
    geno_sim["node"] = get_set(arg)
    founders = get_founders(arg)
    founder_nodes = founder_nodes[sites]
    
    for i in range(0, len(founders)):
        geno_sim.loc[geno_sim["node"] == founders[i], sites] = list(founder_nodes.iloc[founders[i]])

    edges = get_edges(founders, arg)

    while edges.empty != True:
        for i in range(0, len(edges["parent"])):
            snps = genmap.loc[(genmap['Position(bp)'] >= edges.iloc[i]["left"]) & (genmap['Position(bp)'] <= edges.iloc[i]["right"]), "Marker"]
            geno_sim.loc[geno_sim["node"] == edges.iloc[i]["child"], snps] = geno_sim.loc[geno_sim["node"] == edges.iloc[i]["parent"], snps].values[0]
        edges = get_edges(edges["child"], arg)
    return(geno_sim)

#return offspring genotypes of genotype propagation
def get_offspring_geno(arg: tskit.TreeSequence, geno_sim: pd.DataFrame) -> pd.DataFrame:
    """
    Get genotypes of offspring nodes from propagated genotypes DataFrame.

    Parameters:
        arg (tskit.TreeSequence): TreeSequence object obtained from msprime.sim_ancestry.
        geno_sim (pd.DataFrame): DataFrame containing propagated genotypes.

    Returns:
        pd.DataFrame: DataFrame containing genotypes of offspring.
    """
    return(geno_sim.loc[geno_sim["node"].isin(get_offspring(arg))])

#read_hapmap function to read recombination rates
def get_rate_map(genmap: pd.DataFrame) -> msprime.RateMap:
    """
    Get a rate map from a DataFrame.

    Parameters:
        genmap (pd.DataFrame): DataFrame containing rate map information.

    Returns:
        msprime.RateMap: Generated rate map.
    """
    genmap = genmap[["Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)"]]
    return(msprime.RateMap.read_hapmap(io.StringIO(genmap.to_string(index = False)), rate_col = 2))
    
#function to join nodes at individual level to reconstruct diploid genome
def join_nodes(arg: tskit.TreeSequence, geno_sim: pd.DataFrame) -> pd.DataFrame:
    """
    Join nodes at individual level to reconstruct diploid genome.

    Parameters:
        arg (tskit.TreeSequence): TreeSequence object obtained from msprime.sim_ancestry.
        geno_sim (pd.DataFrame): DataFrame containing propagated genotypes.

    Returns:
        pd.DataFrame: DataFrame containing reconstructed diploid genomes.
    """
    #get nodes and individuals from node table
    offspring = get_offspring(arg)
    individual = []
    for i in range(0, len(offspring)):
        individual.append(arg.node(offspring[i]).individual)
    #merge nodes/ind table with geno_sim nodes
    nodes_ind = pd.DataFrame({"node": offspring, "individual": individual})
    offspring_geno = pd.merge(nodes_ind, geno_sim, on = "node", how = "inner")
    #aggregate nodes per individual
    offspring_geno = offspring_geno.drop(["node"], axis = 1).groupby("individual").agg(lambda x: "".join(x)).reset_index()
    return(offspring_geno)
        
#function to recast snp into 0,1,2 depending on encoding (in this case: B73 allele: -1, heterozygote 0, non-B73 allele: 1)
def additive_encoding(ref: pd.DataFrame, genotypes: pd.DataFrame) -> pd.DataFrame:
    """
    Recast SNP into 0, 1, 2 depending on encoding.

    Parameters:
        ref (pd.DataFrame): DataFrame containing reference alleles.
        genotypes (pd.DataFrame): DataFrame containing genotypes.

    Returns:
        pd.DataFrame: DataFrame containing recoded genotypes.
    """
    #melt genotypes, merge to ref allele by snp, for each SNP if heterozygous -> 1 else if identical to b73 -> 0 else 2 
    melted = genotypes.melt(id_vars=["individual"], var_name = "SNP", value_name = "sim_allele")
    merged = pd.merge(ref, melted, on = "SNP", how = "inner")
    merged["Match"] = merged["B73"] == merged["sim_allele"]
    merged["Value"] = np.where(merged["sim_allele"].str[0] == merged["sim_allele"].str[1], np.where(merged["Match"] == True, -1, 1), 0)
    geno_add = merged.pivot(index = "individual", columns = "SNP", values = "Value").reset_index()
    
    return(geno_add)

def genotype_simulation(genetic_map: pd.DataFrame, parent_genos: pd.DataFrame, ref_allele: pd.DataFrame,
                        founder_list: List[str], offspring: int, selfing_genos: int = 5) -> (pd.DataFrame, List[tskit.TreeSequence]):
    """
    Simulate genotypes based on the given parameters.

    Parameters:
        genetic_map (pd.DataFrame): DataFrame containing genetic map information.
        parent_genos (pd.DataFrame): DataFrame containing genotypes of parent individuals.
        ref_allele (pd.DataFrame): DataFrame containing reference alleles.
        founder_list (List[str]): List of two founder names.
        offspring (int): Number of offspring.
        selfing_genos (int): Number of selfing generations.

    Returns:
        pd.DataFrame: DataFrame containing simulated genotypes.
        List[tskit.TreeSequence]: List containing ARGs of the simulated pedigree per chromosome.
    """
    #get founder nodes
    founder_nodes = get_founder_nodes(parent_genos, founder_list)
    #reduce ref_alleles to alleles in genmap
    ref_allele = ref_allele[ref_allele["SNP"].isin(genetic_map["Marker"])]
    #set up pedigree
    pedigree = cross_selfing_ped(offspring = offspring, selfing_genos = selfing_genos)
    #init final genotypes dataframe
    genotypes = pd.DataFrame({"individual": pedigree.loc[pedigree["time"] == 0, "id"]})
    #init final ARG list
    ARGs = list()
    #loop over each chromosome, append results (treats chromosomes entirely independent)
    #possible to parallelise due to independence of tasks
    print("starting simulation")
    for i in genetic_map["Chromosome"].unique():
        #reduce genetic map to chr_i
        chr_genmap = genetic_map.loc[genetic_map["Chromosome"] == i]
        #set up rate map for chr_i
        chr_rate_map = get_rate_map(chr_genmap)
        #get chr_i length
        chr_len = chr_genmap["Position(bp)"].max()
        #turn pedigree into ts with chr_i length
        chr_ts = df_to_ts(pedigree, seq_len=chr_len)
        #simulate chr_i ARG
        chr_arg = msprime.sim_ancestry(initial_state = chr_ts, model="fixed_pedigree", recombination_rate = chr_rate_map)
        #propagate ARG recombinations to offspring genotypes
        chr_geno_sim = propagate_geno(chr_arg, founder_nodes, chr_genmap)
        #join haploid offspring simulation nodes to diploid individuals
        chr_genotypes = join_nodes(chr_arg, chr_geno_sim)
        #merge chr_i genotypes to final genotypes
        genotypes = pd.merge(genotypes, chr_genotypes, on = "individual", how = "inner")
        #append ARG to arg list
        ARGs.append(chr_arg)
        print(f"finished chromsome {i}")

    

    return(genotypes, ARGs)