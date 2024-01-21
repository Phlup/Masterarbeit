from project_lib.genotype_simulation import *
from project_lib.stat_functions import *

# simulates offspring genotypes of randomly chosen NAM parents (no real reference) to use in phenotype prediction model
# usage: Read in genetic map, parental genotypes and reference allele to simulate genotypes of offspring w.r.t. params
# set in genotype_simulation()

if __name__ == '__main__':
    # read genetic map
    genmap = pd.read_csv("data/sim_data/B73_genmap.csv")
    # read parent genotypes
    parent_genos = pd.read_csv("data/sim_data/NAM_parent_genos.csv")
    # read reference allele
    ref_allele = pd.read_csv("data/sim_data/B73_alleles.csv")
    # build artificial pop dict
    populations = pd.read_csv("data/sim_data/populations.csv")
    sim_parents = list()
    for i in range(0, 200):
        cross = populations["parent"].sample(n=2).reset_index(drop=True)
        cross = cross[0] + "_" + cross[1]
        sim_parents.append(cross)
    sim_parents = list(set(sim_parents))
    sim_populations = pd.DataFrame({"pop": range(1, len(sim_parents) + 1), "name": sim_parents})

    sim_pop_dict = dict(zip(sim_populations["pop"], sim_populations["name"]))
    sim_pops = sim_populations["pop"]

    for i in sim_pops:
        pop = sim_pop_dict[i].split("_")

        # normal estimated recombination rate
        norm_sim = genotype_simulation(genetic_map=genmap, parent_genos=parent_genos, ref_allele=ref_allele,
                                       founder_list=pop, offspring=200, selfing_genos=5)

        norm_sim[0].to_csv("sim_output/parent_sim/geno_encoding/geno_" + str(i) + ".csv", index=False)

        norm_add = additive_encoding(ref_allele, norm_sim[0])

        norm_add.to_csv("sim_output/parent_sim/additive_encoding/add_" + str(i) + ".csv", index=False)
        print(f"finished simulating {sim_pop_dict[i]} ({i}/{len(sim_pop_dict)})")
    sim_populations.to_csv("sim_output/parent_sim/sim_populations.csv", index = False)