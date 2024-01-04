from project_lib.genotype_simulation import *

#usage: Read in genetic map, parental genotypes and reference allele to simulate genotypes of offspring w.r.t. params
#set in genotype_simulation()
#simulates all reference pop offspring genotypes and saves ref pop as additive encoding

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
    pops = populations["pop"]

    #genmap with high recombination rate
    genmap_high = genmap.copy()
    genmap_high["Rate(cM/Mb)"] += 5
    #genmap with 0 recombination rate
    genmap_zero = genmap.copy()
    genmap_zero["Rate(cM/Mb)"] = 0
    #genmap with constant recombination rate (mean over all recombination)
    genmap_mean = genmap.copy()
    genmap_mean["Rate(cM/Mb)"] = genmap["Rate(cM/Mb)"].mean()

    recomb_stats = pd.DataFrame({"pop": range(0, 12), "mean_norm": None, "num_norm": None, "mean_high": None, "num_high": None,
         "mean_zero": None, "num_zero": None, "mean_mean": None, "num_mean": None})

    for i in pops:
        pop = pop_dict[i].split("_")

        #normal estimated recombination rate
        norm_sim = genotype_simulation(genetic_map=genmap, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        norm_sim.to_csv("sim_output/normal_rec/geno_encoding/geno_" + str(i) + ".csv", index = False)

        norm_add = additive_encoding(ref_allele, norm_sim)

        norm_add.to_csv("sim_output/normal_rec/additive_encoding/add_" + str(i) + ".csv", index = False)

        #high recombination rate
        high_sim = genotype_simulation(genetic_map=genmap_high, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        high_sim.to_csv("sim_output/high_rec/geno_encoding/geno_" + str(i) + ".csv", index = False)

        high_add = additive_encoding(ref_allele, high_sim)

        high_add.to_csv("sim_output/high_rec/additive_encoding/add_" + str(i) + ".csv", index = False)

        #zero recombination rate
        zero_sim = genotype_simulation(genetic_map=genmap_zero, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        zero_sim.to_csv("sim_output/zero_rec/geno_encoding/geno_" + str(i) + ".csv", index = False)

        zero_add = additive_encoding(ref_allele, zero_sim)

        zero_add.to_csv("sim_output/zero_rec/additive_encoding/add_" + str(i) + ".csv", index = False)

        #mean recombination rate
        mean_sim = genotype_simulation(genetic_map=genmap_mean, parent_genos=parent_genos, ref_allele=ref_allele,
                                        founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        mean_sim.to_csv("sim_output/mean_rec/geno_encoding/geno_" + str(i) + ".csv", index = False)

        mean_add = additive_encoding(ref_allele, mean_sim)

        mean_add.to_csv("sim_output/mean_rec/additive_encoding/add_" + str(i) + ".csv", index = False)
        print(f"finished simulating {pop_dict[i]} ({i}/{len(pop_dict)})")

    #save all real pop genos as add encoding
    for i in pops:
        real_i = pd.read_csv("data/NAM_genotype_data/geno_encoding/pop_" + str(i) + "_genos.csv")
        add_i = additive_encoding(ref_allele, real_i)
        add_i.to_csv("data/NAM_genotype_data/additive_encoding/pop_" + str(i) + "_add.csv", index = False)