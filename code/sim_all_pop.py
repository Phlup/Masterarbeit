from project_lib.genotype_simulation import *
from project_lib.stat_functions import *

# simulates all reference pop offspring genotypes, calculates recombination stats and saves ref pop as additive encoding
# usage: Read in genetic map, parental genotypes and reference allele to simulate genotypes of offspring w.r.t. params
# set in genotype_simulation()

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

    genmap_high = genmap.copy()
    genmap_high["Rate(cM/Mb)"] += 5
    # ensure last recombination rate is 0 at every chromosome (requirement for msprime)
    genmap_high.loc[genmap_high.groupby("Chromosome").tail(1).index, "Rate(cM/Mb)"] = 0
    # genmap with 0 recombination rate
    genmap_zero = genmap.copy()
    genmap_zero["Rate(cM/Mb)"] = 0
    # genmap with constant recombination rate (mean over all recombination)
    genmap_mean = genmap.copy()
    genmap_mean["Rate(cM/Mb)"] = genmap["Rate(cM/Mb)"].mean()
    genmap_mean.loc[genmap_mean.groupby("Chromosome").tail(1).index, "Rate(cM/Mb)"] = 0

    recomb_stats = pd.DataFrame({"pop": pops, "mean_norm": None, "num_norm": None, "mean_high": None, "num_high": None,
                                 "mean_zero": None, "num_zero": None, "mean_mean": None, "num_mean": None})

    for i in pops:
        pop = pop_dict[i].split("_")

        # normal estimated recombination rate
        norm_sim = genotype_simulation(genetic_map=genmap, parent_genos=parent_genos, ref_allele=ref_allele,
                                       founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        recomb_stat = calc_total_recomb(norm_sim[1])
        recomb_stats.loc[recomb_stats["pop"] == i, "mean_norm"] = recomb_stat[0]
        recomb_stats.loc[recomb_stats["pop"] == i, "num_norm"] = recomb_stat[1]

        norm_sim[0].to_csv("sim_output/normal_rec/geno_encoding/geno_" + str(i) + ".csv", index=False)

        norm_add = additive_encoding(ref_allele, norm_sim[0])

        norm_add.to_csv("sim_output/normal_rec/additive_encoding/add_" + str(i) + ".csv", index=False)

        # high recombination rate
        high_sim = genotype_simulation(genetic_map=genmap_high, parent_genos=parent_genos, ref_allele=ref_allele,
                                       founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        recomb_stat = calc_total_recomb(high_sim[1])
        recomb_stats.loc[recomb_stats["pop"] == i, "mean_high"] = recomb_stat[0]
        recomb_stats.loc[recomb_stats["pop"] == i, "num_high"] = recomb_stat[1]

        high_sim[0].to_csv("sim_output/high_rec/geno_encoding/geno_" + str(i) + ".csv", index=False)

        high_add = additive_encoding(ref_allele, high_sim[0])

        high_add.to_csv("sim_output/high_rec/additive_encoding/add_" + str(i) + ".csv", index=False)

        # zero recombination rate
        zero_sim = genotype_simulation(genetic_map=genmap_zero, parent_genos=parent_genos, ref_allele=ref_allele,
                                       founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        recomb_stat = calc_total_recomb(zero_sim[1])
        recomb_stats.loc[recomb_stats["pop"] == i, "mean_zero"] = recomb_stat[0]
        recomb_stats.loc[recomb_stats["pop"] == i, "num_zero"] = recomb_stat[1]

        zero_sim[0].to_csv("sim_output/zero_rec/geno_encoding/geno_" + str(i) + ".csv", index=False)

        zero_add = additive_encoding(ref_allele, zero_sim[0])

        zero_add.to_csv("sim_output/zero_rec/additive_encoding/add_" + str(i) + ".csv", index=False)

        # mean recombination rate
        mean_sim = genotype_simulation(genetic_map=genmap_mean, parent_genos=parent_genos, ref_allele=ref_allele,
                                       founder_list=pop, offspring=size_dict[i], selfing_genos=5)

        recomb_stat = calc_total_recomb(mean_sim[1])
        recomb_stats.loc[recomb_stats["pop"] == i, "mean_mean"] = recomb_stat[0]
        recomb_stats.loc[recomb_stats["pop"] == i, "num_mean"] = recomb_stat[1]

        mean_sim[0].to_csv("sim_output/mean_rec/geno_encoding/geno_" + str(i) + ".csv", index=False)

        mean_add = additive_encoding(ref_allele, mean_sim[0])

        mean_add.to_csv("sim_output/mean_rec/additive_encoding/add_" + str(i) + ".csv", index=False)
        print(f"finished simulating {pop_dict[i]} ({i}/{len(pop_dict)})")

    # save recombination stats
    recomb_stats.to_csv("stats/summary_stats/recomb_stats.csv", index=False)
    # save all real pop genos as add encoding
    for i in pops:
        real_i = pd.read_csv("data/NAM_genotype_data/geno_encoding/pop_" + str(i) + "_genos.csv")
        add_i = additive_encoding(ref_allele, real_i)
        add_i.to_csv("data/NAM_genotype_data/additive_encoding/pop_" + str(i) + "_add.csv", index=False)
