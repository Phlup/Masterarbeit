from project_lib.stat_functions import *


if __name__ == '__main__':
    # read genetic map
    genmap = pd.read_csv("data/sim_data/B73_genmap.csv")
    # read population names
    populations = pd.read_csv("data/sim_data/populations.csv")
    pops = populations["pop"]
    #vary over recombination scenarios
    rec_param = list(["normal_rec", "high_rec", "zero_rec", "mean_rec"])

    for i in rec_param:
        for j in pops:
            sim_add = pd.read_csv("sim_output/" + i + "/additive_encoding/add_" + str(j) + ".csv")
            real_add = pd.read_csv("data/NAM_genotype_data/additive_encoding/pop_" + str(j) + "_add.csv")
            sim_add = sim_add[sim_add.columns.intersection(genmap["Marker"])]
            real_add = real_add[real_add.columns.intersection(genmap["Marker"])]
            summary_plot(real_additive=real_add, sim_additive=sim_add,
                         pop_num=j, out_path="plots/dist_plots/" + i + "/summary_" + str(j) + ".png")