import tskit
import msprime
import stdpopsim
import networkx as nx
import matplotlib.pyplot as plt
import io
from IPython.display import display as SVG
import sys
import tszip
import pandas as pd
import numpy as np

#draw pedigree
def draw_pedigree(ped_ts):
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
    plt.show()

#print ancestry
def draw_ancestry(ts):
    node_labels = {node.id: f"{node.individual}({node.id})" for node in ts.nodes()}
    SVG(ts.draw_svg(y_axis=True,  node_labels=node_labels, size=(3500,400)))

#generate pedigree df with n "founder" individuals
def pedigree_init(n = 10):
    pedigree = pd.DataFrame(
    {"id": [i for i in range(0, n)],
    "parent0": ["."]*n,
    "parent1": ["."]*n,
    "time": [0]*n})
    return pedigree

#generate msprime trees from df
def df_to_ts(df, seq_len = 100):
    ped_string = df.to_string(index = False)
    ts_ped = msprime.parse_pedigree(io.StringIO("#" + ped_string), sequence_length = seq_len)
    return ts_ped

#selfing (adding selfed individual to pedigree)
def add_selfing(df, size_diff = 0):
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
def add_random_mating(df, size_diff = 0):
    df['time'] += 1
    parents = df.loc[df['time'] == df['time'].min(), "id"].to_numpy()
    max_id = df.loc[df['id'].max(), "id"]
    ids = [i for i in range(max_id + 1, max_id + 1 + len(parents) + size_diff)]
    df = pd.concat(objs=(df, pd.DataFrame({"id" : ids, "parent0": np.random.choice(parents, len(parents) + size_diff),
                                           "parent1": np.random.choice(parents, len(parents) + size_diff),
                                           "time" : [0]*(len(parents) + size_diff)}))).reset_index(drop = True)
    
    return(df)
    
#selective mating
def add_selective_mating(df, parents, offspring = 1):
    df['time'] += 1
    max_id = df.loc[df['id'].max(), "id"]
    ids = [i for i in range(max_id + 1, max_id + 1 + offspring)]
    df = pd.concat(objs=(df, pd.DataFrame({"id": ids, "parent0": np.resize(parents[0], offspring),
                                           "parent1": np.resize(parents[1], offspring),
                                           "time": [0]*offspring}))).reset_index(drop = True)     
    return(df)


#propagate geno functions
def get_set(ts):
    return(list(set(ts.edges_parent).union(set(ts.edges_child))))

def get_founders(ts):
    return(list(set(ts.edges_parent) - set(ts.edges_child)))

def get_offspring(ts):
    return(list(set(ts.edges_child) - set(ts.edges_parent)))

def get_edges(parents, ts):
    pc_df = pd.DataFrame({"left": ts.edges_left, "right": ts.edges_right, "parent" : ts.edges_parent, "child" : ts.edges_child})
    children = pc_df.loc[pc_df["parent"].isin(parents),]
    return(children)

#replace snp_ids with df of snp + position
def propagate_geno(ts, genotypes, snp_ids):
    geno_sim = pd.DataFrame(columns = list(["node"]) + list(snp_ids))
    geno_sim["node"] = get_set(ts)
    founders = get_founders(ts)
    sites = len(pure_genos.columns)
    for i in range(0, len(founders)):
        geno_sim.loc[geno_sim["node"] == founders[i], pure_genos.columns] = list(pure_genos.iloc[founders[i]])

    edges = get_edges(founders,ts)

    while edges.empty != True:
        for i in range(0, len(edges["parent"])):
            #replace SNPS with function that takes left and right and SNP names between them from SNP+position df
            snps = pure_genos.columns[round((edges.iloc[i]["left"]/100)*sites):round((edges.iloc[i]["right"]/100)*sites)]
            geno_sim.loc[geno_sim["node"] == edges.iloc[i]["child"], snps] = geno_sim.loc[geno_sim["node"] == edges.iloc[i]["parent"], snps].values[0]
        edges = get_edges(edges["child"], ts)
    
    
    return(geno_sim)
