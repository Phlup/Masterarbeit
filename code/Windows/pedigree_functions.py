import tskit
import msprime
import stdpopsim
import networkx as nx
import matplotlib.pyplot as plt
import io
from IPython.display import display as SVG
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

#function to choose founders and split diploid genome into ts nodes -> output genoypes df to use in propagate_geno function
#not generalizable, specific to dataset from bukler et al 2007
def get_founder_nodes(genotypes, founders):
    founder_geno = genotypes.loc[genotypes['RIL'].isin(founders)].drop("RIL", axis = 1)
    founder_nodes = pd.concat([founder_geno.loc[0].str[0],
                              founder_geno.loc[0].str[2],
                              founder_geno.loc[1].str[0],
                              founder_geno.loc[1].str[2]], axis = 1).T.reset_index(drop = True)
    return founder_nodes

#add typing!!
def propagate_geno(arg, founder_nodes, genmap):
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

def get_offspring_geno(arg, geno_sim):
    return(geno_sim.loc[geno_sim["node"].isin(get_offspring(arg))])

#read_hapmap function to read recombination rates
def get_rate_map(genmap):
    genmap = genmap[["Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)"]]
    return(msprime.RateMap.read_hapmap(io.StringIO(genmap.to_string(index = False))))
    
#function to join nodes at individual level to reconstruct diploid genome
def join_nodes(arg, geno_sim):
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
    
#function to recast snp into 0,1,2 depending on encoding (in this case: B73 allele: 0, heterozygote 1, non-B73 allele: 2)
def additive_encoding(ref, genotypes):
    #melt genotypes, merge to ref allele by snp, for each SNP if heterozygous -> 1 else if identical to b73 -> 0 else 2 
    melted = genotypes.melt(id_vars=["individual"], var_name = "SNP", value_name = "sim_allele")
    merged = pd.merge(ref, melted, on = "SNP", how = "inner")
    merged["Match"] = merged["B73"] == merged["sim_allele"]
    merged["Value"] = np.where(merged["sim_allele"].str[0] == merged["sim_allele"].str[1], np.where(merged["Match"] == True, 0, 2), 1)
    geno_add = merged.pivot(index = "individual", columns = "SNP", values = "Value")
    
    return(geno_add)