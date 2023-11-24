from project_lib.genotype_simulation import *

#usage: Read in genetic map, parental genotypes and reference allele to simulate genotypes of offspring w.r.t. params
#set in genotype_simulation()

if __name__ == '__main__':
    # read genetic map
    genmap = pd.read_csv("../data/test_data/B73_genmap.csv")
    # read parent genotypes
    parent_genos = pd.read_csv("../data/test_data/NAM_parent_genos.csv")
    # read reference allele
    ref_allele = pd.read_csv("../data/test_data/B73_alleles.csv")

    whole_sim = genotype_simulation(genetic_map=genmap, parent_genos=parent_genos, ref_allele=ref_allele,
                                    founder_list=list(["B73", "B97"]), offspring=194, selfing_genos=5)

    whole_sim.to_csv("../sim_output/whole_sim.csv")