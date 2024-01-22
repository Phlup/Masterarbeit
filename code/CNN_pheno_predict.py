import pandas as pd
import tensorflow as tf
from tensorflow.keras import models, layers
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# trains CNN on parental traits + genetic map features of simulated populations to predict trait summary stats
# on real populations, scores RMSE, correlation(predictions, test_y) and if best performing population was accurately predicted

if __name__ == '__main__':
    #load train and test data (marker traits/genmap features/correlations on sim and real populations)
    real_parent_trait = pd.read_csv("stats/pheno_prediction/real_parent_trait_map.csv")
    sim_parent_trait = pd.read_csv("stats/pheno_prediction/sim_parent_trait_map.csv")

    real_parent_cors = pd.read_csv("stats/pheno_prediction/real_cors.csv")
    sim_parent_cors = pd.read_csv("stats/pheno_prediction/sim_cors.csv")

    real_phenos = pd.read_csv("stats/pheno_prediction/real_summary.csv")
    sim_phenos = pd.read_csv("stats/pheno_prediction/sim_pops_summary.csv")

    #load genmap posdiff data
    genmap_posdiff = pd.read_csv("data/sim_data/genmap_posdiff.csv")
    mean_posdiff = genmap_posdiff["pos_diff"].mean() / 1000000
    #load recombination sumstats to calc window size
    geno_sumstats = pd.read_csv("stats/summary_stats/sum_results.csv")
    mean_recomb_length = geno_sumstats[geno_sumstats["rec_param"] == "normal_rec"]["recomb_mean"]
    #kernel_y in CNN as mean markers per recombination block
    window_size = int(round(mean_recomb_length/mean_posdiff).iloc[0])
    print(window_size)
    #set features for CNN (out of p1, p2, position, rate, cM, pos_diff)
    feature_combinations = list()
    feature_combinations.append(["p1", "p2"])
    feature_combinations.append(["p1", "p2", "rate"])
    feature_combinations.append(["p1", "p2", "cM"])
    feature_combinations.append(["p1", "p2", "position"])
    feature_combinations.append(["p1", "p2", "pos_diff"])
    feature_combinations.append(["p1", "p2", "pos_diff","rate"])
    feature_combinations.append(["corr"])
    ##output df
    pred_results = pd.DataFrame(columns=["model", "trait", "task", "rmse", "corr", "corr_p", "best_predict"])
    for features in feature_combinations:

        real_parent_trait_features = real_parent_trait[real_parent_trait["feature"].isin(features)].drop(columns=["feature"])
        sim_parent_trait_features = sim_parent_trait[sim_parent_trait["feature"].isin(features)].drop(columns=["feature"])
        input_y = len(genmap_posdiff)
        if (features[0]) == "corr":
            real_parent_trait_features = real_parent_cors
            sim_parent_trait_features = sim_parent_cors
            input_y = real_parent_cors.drop(columns=["trait"]).shape[1]

        ##describe CNN model
        model = models.Sequential()
        #input_dims = (#num_rows_per_sample, #num_columns_per_sample)
        #kernel = (#num_trait_map_features, #window_size)
        model.add(layers.Conv2D(32, (len(features), window_size), activation='relu', input_shape=(len(features), input_y, 1)))
        model.add(layers.MaxPooling2D((1, 4)))
        model.add(layers.Flatten())
        model.add(layers.Dense(32, activation='relu'))
        model.add(layers.Dense(1, activation='linear'))

        #compile model
        model.compile(optimizer='adam', loss='mean_squared_error')
        traits = list(["silk", "tassel", "oil", "protein", "starch"])
        tasks = list(["trait_mean", "trait_95_perc"])
        for i in traits:
            for j in tasks:
                real_x = real_parent_trait_features[real_parent_trait_features["trait"].isin([i])].drop(columns=["trait"])
                real_x = np.array(real_x).reshape((int(len(real_x)/len(features)), len(features), input_y, 1))
                sim_x = sim_parent_trait_features[sim_parent_trait_features["trait"].isin([i])].drop(columns=["trait"])
                sim_x = np.array(sim_x).reshape((int(len(sim_x)/len(features)), len(features), input_y, 1))
                real_y = real_phenos[real_phenos["trait"].isin([i])][j]
                sim_y = sim_phenos[sim_phenos["trait"].isin([i])][j]
                model.fit(sim_x, sim_y, epochs=500, batch_size=1)
                predictions = model.predict(real_x)
                loss = model.evaluate(real_x, real_y)
                rmse = np.sqrt(loss)
                correlation, corr_p = pearsonr(predictions.flatten(), real_y)
                best_pred = (np.argmax(predictions.flatten()) == np.argmax(real_y))
                pred_results = pred_results._append({"model": '_'.join(features), "trait": i, "task": j, "rmse": round(rmse,3),
                                              "corr": round(correlation,3),
                                              "corr_p": round(corr_p,3 ) if corr_p >= 0.001 else "<0.001", "best_predict": best_pred},
                                             ignore_index=True)
                plt.scatter(predictions.flatten(), real_y)
                plt.xlabel("Predicted phenotype")
                plt.ylabel("True phenotype")
                plt.title("")
                plt.savefig("plots/prediction_plots/CNN/CNN_"+i+"_"+j+"_"+'_'.join(features)+".png", dpi = 500, bbox_inches = "tight")
                plt.close()
    pred_results.to_csv("stats/pheno_prediction/results/pred_results_CNN.csv", index = False)


