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

    founders = parent_genos["RIL"]

    founder_pairs = [(founders[0], non_b73) for non_b73 in founders[1:]]

    for i in range(0, len(founder_pairs)):

        whole_sim = genotype_simulation(genetic_map=genmap, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=founder_pairs[i], offspring=200, selfing_genos=5)

        whole_sim.to_csv("../sim_output/geno_encoding/geno_" + founder_pairs[i][0] + "_" + founder_pairs[i][1] + ".csv")

        add_sim = additive_encoding(ref_allele, whole_sim)

        add_sim.to_csv("../sim_output/additive_encoding/geno_" + founder_pairs[i][0] + "_" + founder_pairs[i][1] + ".csv")
        print("finished simulating " + founder_pairs[i][0] + "x" + founder_pairs[i][1] + " (" + str(i+1) + "/"
              + str(len(founder_pairs)) + ")")

