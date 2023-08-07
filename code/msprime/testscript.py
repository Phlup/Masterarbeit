import tskit
import msprime
import stdpopsim

# import sunflower species
from svgwrite.container import SVG

species = stdpopsim.get_species("HelAnn")

# create selfing pedigree


pb = msprime.PedigreeBuilder()
mom_id = pb.add_individual(time=5)
dad_id = pb.add_individual(time=5)

for counter in range(0, 5):
    # as this are diploids and the crops are all homozygous we mimic this using one mom instead of a mom and a dad
    # but this would equal creating a clone of the mother and mating the mother and her clone
    # this would create less heterozygotes as expected as allele 1 of the mother could recombine with allele 1 of the clone and they are identical thus not generating any heterozygotes
    # if the father is a different homozygous sequence than each recombination should also create heterozygotes
    # we just have to postprocess the sequence as the alleles 0 and 1 of the mom_id and the alleles 2 and 3 of the dad_id are identical
    # maybe I can add an artificial merger event using msprime?
    previous_id = pb.add_individual(time=4, parents=[mom_id, dad_id])
    for gen in reversed(range(0, 4)):
        pb.add_individual(time=gen, parents=[previous_id, previous_id], is_sample=(gen == 0))


pedigree = pb.finalise(sequence_length=100)
# TODO replace with display(pedigree) when its implemented in tskit
# https://github.com/tskit-dev/tskit/issues/2093
print(pedigree)


ped_ts = msprime.sim_ancestry(
    initial_state=pedigree, model="fixed_pedigree", random_seed=42,
    recombination_rate=0.001)
b
print(ped_ts.draw_text())






# IBD Segments
segments = ped_ts.ibd_segments(store_pairs=True)
print(list(segments.keys()))

for pair, value in segments.items():
    print(pair, "::", value)



segments = ped_ts.ibd_segments(store_pairs=True, store_segments=True)
#segments[(0, 1)]

import networkx as nx
import matplotlib.pyplot as plt


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


# this is an alternative way to read in pedigrees
import io

ped_txt = """\
# id parent0 parent1 time is_sample
0   5   4   0.0 1
1   3   3   0.0 1
2   3   4   0.0 1
3   7   7   1.0 0
4   8   8   1.0 0
5   8   6   1.0 0
6   9   10  2.0 0
7   10  10  2.0 0
8   9   9   2.0 0
9   11  12  3.0 0
10  11  12  3.0 0
11  .   .   4.0 0
12  .   .   4.0 0
"""
pedigree2 = msprime.parse_pedigree(io.StringIO(ped_txt), sequence_length=100)

draw_pedigree(pedigree2.tree_sequence())







ped_ts = msprime.sim_ancestry(
    initial_state=pedigree, model="fixed_pedigree", random_seed=41,
    recombination_rate=0.001)
node_labels = {node.id: f"{node.individual}({node.id})" for node in ped_ts.nodes()}
#SVG(ped_ts.draw_svg(y_axis=True, node_labels=node_labels, size=(600, 200)))