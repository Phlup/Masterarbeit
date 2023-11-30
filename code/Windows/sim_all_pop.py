from project_lib.genotype_simulation import *

#usage: Read in genetic map, parental genotypes and reference allele to simulate genotypes of offspring w.r.t. params
#set in genotype_simulation()

if __name__ == '__main__':
    # read genetic map
    genmap = pd.read_csv("data/sim_data/B73_genmap.csv")
    # read parent genotypes
    parent_genos = pd.read_csv("data/sim_data/NAM_parent_genos.csv")
    # read reference allele
    ref_allele = pd.read_csv("data/sim_data/B73_alleles.csv")
    # read population names
    populations = pd.read_csv("data/sim_data/populations.csv")

    pop_dict = dict(zip(populations["pop"], populations["name"]))
    size_dict = dict(zip(populations["pop"], populations["size"]))

    #genmap with high recombination rate
    genmap_high = genmap.copy()
    genmap_high["Rate(cM/Mb)"] += 5
    #genmap with 0 recombination rate
    genmap_zero = genmap.copy()
    genmap_zero["Rate(cM/Mb)"] = 0
    #genmap with constant recombination rate (mean over all recombination)
    genmap_mean = genmap.copy()
    genmap_mean["Rate(cM/Mb)"] = genmap["Rate(cM/Mb)"].mean()

    for i in range(0, len(pop_dict)):
        pop = pop_dict[i+1].split("_")

        #normal estimated recombination rate
        norm_sim = genotype_simulation(genetic_map=genmap, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i+1], selfing_genos=5)

        norm_sim.to_csv("sim_output/normal_rec/geno_encoding/geno_" + str(i+1) + ".csv")

        norm_add = additive_encoding(ref_allele, norm_sim)

        norm_add.to_csv("sim_output/normal_rec/additive_encoding/add_" + str(i+1) + ".csv")

        #high recombination rate
        high_sim = genotype_simulation(genetic_map=genmap_high, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i+1], selfing_genos=5)

        high_sim.to_csv("sim_output/high_rec/geno_encoding/high_geno_" + str(i+1) + ".csv")

        high_add = additive_encoding(ref_allele, high_sim)

        high_add.to_csv("sim_output/high_rec/additive_encoding/high_add_" + str(i+1) + ".csv")

        #zero recombination rate
        zero_sim = genotype_simulation(genetic_map=genmap_zero, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i+1], selfing_genos=5)

        zero_sim.to_csv("sim_output/zero_rec/geno_encoding/zero_geno_" + str(i+1) + ".csv")

        zero_add = additive_encoding(ref_allele, zero_sim)

        zero_add.to_csv("sim_output/zero_rec/additive_encoding/zero_add_" + str(i+1) + ".csv")

        #mean recombination rate
        mean_sim = genotype_simulation(genetic_map=genmap_mean, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i+1], selfing_genos=5)

        mean_sim.to_csv("sim_output/mean_rec/geno_encoding/mean_geno_" + str(i+1) + ".csv")

        mean_add = additive_encoding(ref_allele, mean_sim)

        mean_add.to_csv("sim_output/mean_rec/additive_encoding/mean_add_" + str(i+1) + ".csv")
        print(f"finished simulating {pop_dict[i+1]} ({i+1}/{len(pop_dict)})")

