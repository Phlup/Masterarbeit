import tskit
import msprime
import networkx as nx
import matplotlib.pyplot as plt
import io
from IPython.display import display as SVG
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict
import scipy.stats as stats
import seaborn as sns

#draw pedigree
def draw_pedigree(ped_ts: tskit.TableCollection, out_path = str) -> None:
    """
    Draw a pedigree using NetworkX and matplotlib.

    Parameters:
        ped_ts (tskit.TableCollection): A table collection representing the pedigree.
        out_path (str): Output path for the pedigree plot.

    Returns:
        None, displays plot
    """
    G = nx.DiGraph()
    for ind in ped_ts.individuals():
        time = ped_ts.node(ind.nodes[0]).time
        pop = ped_ts.node(ind.nodes[0]).population
        G.add_node(ind.id, time=time, population=pop)
        for p in ind.parents:
            if p != tskit.NULL:
                G.add_edge(ind.id, p)
    pos = nx.multipartite_layout(G, subset_key="time", align="horizontal")
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    node_colours = [colours[node_attr["population"]] for node_attr in G.nodes.values()]
    nx.draw_networkx(G, pos, with_labels=True, node_color=node_colours)
    plt.savefig(out_path, format="png")

#print ancestry
def draw_ancestry(ts: tskit.TreeSequence, x_size: int = 3000, y_size: int = 400, out_path = None) -> None:
    """
    Print the ancestry of individuals in the tree sequence as ancestral recombination graph (ARG).

    Parameters:
        ts (tskit.TreeSequence): A tree sequence.
        x_size (int): Size of x dimension.
        y_size (int): Size of y dimension.
        out_path (None/str): Output path for the ancestry plot.
    Returns:
        None, displays plot
    """
    node_labels = {node.id: f"{node.individual}({node.id})" for node in ts.nodes()}
    SVG(ts.draw_svg(path = out_path, y_axis=True,  node_labels=node_labels, size=(x_size,y_size)))
    

#calculate recombination tract means and events per ARG
def calc_mean_recomb(arg: tskit.TreeSequence) -> tuple:
    """
    Calculates mean recombination tract length and mean number of recombination events per ARG
    in simulated pedigree offspring 

    Parameters:
        arg (tskit.TreeSequence): TreeSequence object obtained from msprime.sim_ancestry.

    Returns:
        tuple: mean recombination tract length and mean number of recombination events per ARG
    """
    recomb_tracts = pd.DataFrame({"length": abs(arg.edges_left-arg.edges_right), "child": arg.edges_child})
    offspring = list(set(arg.edges_child) - set(arg.edges_parent))
    recomb_tracts = recomb_tracts.loc[recomb_tracts["child"].isin(offspring),]
    return recomb_tracts["length"].mean()/1000000, recomb_tracts["length"].count()/(len(offspring)/2)

#calc total recomb means and events across all chromosomes
def calc_total_recomb(args: List[tskit.TreeSequence]) -> tuple:
    """
    Calculates total recomb means and events across all chromosomes (ARGs) in simulated genotype

    Parameters:
        arg (tskit.TreeSequence): TreeSequence object obtained from msprime.sim_ancestry.

    Returns:
        tuple: total mean recombination tract length and mean number of recombination events per individual (all chromosomes)
    """
    recomb_means = list()
    recomb_nums = list()
    for i in range(0, len(args)):
        recomb_means.append(calc_mean_recomb(args[i])[0])
        recomb_nums.append(calc_mean_recomb(args[i])[1])
    return pd.Series(recomb_means).mean(), pd.Series(recomb_nums).sum()

#calculate allele frequencies
def calc_allele_freq(matrix: List[List[int]], alleles: int = 3) -> List[float]:
    """
    Calculate allele frequencies from a genotype matrix.

    Parameters:
        matrix (List[List[int]]): Genotype matrix.
        alleles (int): Number of alleles per SNP.

    Returns:
        List[float]: List of allele frequencies.
    """
    num_individuals = len(matrix)
    num_snps = len(matrix[0])
    num_alleles = alleles
    allele_frequencies = []

    for j in range(num_snps):
        allele_counts = [0] * num_alleles

        for i in range(num_individuals):
            allele_counts[matrix[i][j]] += 1

        snp_frequencies = [count / num_individuals for count in allele_counts]
        allele_frequencies.extend(snp_frequencies)

    return allele_frequencies

#neis distance
def neis_distance(X: np.ndarray) -> np.ndarray:
    """
    Calculate Nei's genetic distance.

    Parameters:
        X (np.ndarray): Genotype matrix.

    Returns:
        np.ndarray: Nei's genetic distance matrix.
    """
    #matrix product of X and its transpose
    d = np.matmul(X, X.T)
    
    #sqrt of diagonal elements
    vec = np.sqrt(np.diag(d))
    
    #Normalize columns of d
    d /= vec[:, np.newaxis]
    
    #Normalize rows of d
    d /= vec
    
    #negative logarithm
    d = -np.log(d)
    
    #todistance matrix
    d = np.asarray(d)
    
    return d

#tabulate genotype counts across all individuals and sites
def genotype_counts(genotype_sim: pd.DataFrame) -> pd.Series:
    """
    Tabulate genotype counts across all individuals and sites.

    Parameters:
        genotype_sim (pd.DataFrame): DataFrame containing simulated genotypes.

    Returns:
        pd.Series: Series containing genotype counts.
    """
    return(pd.Series(genotype_sim.values.flatten().tolist()).value_counts())

#generate and save summary histogram, KDEs and statistics to compare real and simulated genotype distribution
def summary_plot(real_additive: pd.DataFrame, sim_additive: pd.DataFrame, pop_num: int, out_path: str) -> None:
    """
    Generate a summary plot comparing real and simulated genotypes.

    Parameters:
        real_additive (pd.DataFrame): DataFrame containing real additive encoded genotypes (no individual column!).
        sim_additive (pd.DataFrame): DataFrame containing simulated additive encoded genotypes.
        pop_num (int): number of simulated population
        out_path (str): Output path for the summary plot.
    """
    #get vector of sums of additive encoding for all individuals
    real_sum = real_additive.sum(axis = 1)
    sim_sum = sim_additive.sum(axis = 1)

    #calculate ks test statistic and p value
    ks_statistic, p_value = stats.ks_2samp(real_sum, sim_sum)
    
    #calculate wasserstein distance
    wasserstein_dist = stats.wasserstein_distance(real_sum, sim_sum)

    #compute common bins
    bins = np.linspace(min(np.concatenate([real_sum, sim_sum])), max(np.concatenate([real_sum, sim_sum])), 30)

    #plot histogram and KDEs with ks test and WS dist
    sns.histplot(real_sum, kde = True, label = "Real genotypes", color = "blue", bins = bins, alpha = 0.5)
    sns.histplot(sim_sum, kde = True, label = "Simulated genotypes", color = "orange", bins = bins, alpha = 0.5)
    plt.legend(loc = 'upper right')
    #plt.title("Histograms, KDEs and summary statistics of real and simulated genotypes for population " + str(pop_num))
    plt.xlabel("Sum of additive encoding per individual")
    ks_text = f'KS Test P-Value: {p_value:.3f}' if p_value >= 0.001 else f'KS Test P-Value: <0.001'
    plt.text(0.05, 0.95, ks_text, 
         transform=plt.gca().transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    plt.text(0.05, 0.88, f'Wasserstein Distance: {wasserstein_dist:.2f}', transform=plt.gca().transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.savefig(out_path, dpi = 500, bbox_inches = "tight")
    plt.close()

