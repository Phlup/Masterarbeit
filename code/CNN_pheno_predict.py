import pandas as pd
import tensorflow as tf
from tensorflow.keras import models, layers
import numpy as np
import matplotlib.pyplot as plt

##trains CNN on parental traits + genetic map features of simulated populations to predict trait summary stats
##on real populations, scores RMSE, correlation(predictions, test_y) and if best performing population was accurately predicted

if __name__ == '__main__':
    # load train and test data (marker traits/correlations on sim and real populations)
    real_cors = pd.read_csv("stats/pheno_prediction/real_cors.csv")
    sim_cors = pd.read_csv("stats/pheno_prediction/sim_cors.csv")

    real_parent_trait = pd.read_csv("stats/pheno_prediction/real_parent_trait_map.csv")
    sim_parent_trait = pd.read_csv("stats/pheno_prediction/sim_parent_trait_map.csv")

    real_phenos = pd.read_csv("stats/pheno_prediction/real_summary.csv")
    sim_phenos = pd.read_csv("stats/pheno_prediction/sim_summary.csv")

    # load genmap for marker data
    genmap = pd.read_csv("data/sim_data/B73_genmap.csv")
    ##describe CNN model
    model = models.Sequential()
    # Convolute on (3, 5) window
    input_dims = (#num_rows_per_sample, #num_columns_per_sample)
    kernel = (#num_trait_map_features, #window_size))
    model.add(layers.Conv2D(32, (3, 5), activation='relu', input_shape=(kernel_x, kernel_y, 1)))
    model.add(layers.MaxPooling2D((1, 4)))
    model.add(layers.Flatten())
    model.add(layers.Dense(32, activation='relu'))
    model.add(layers.Dense(1, activation='linear'))

    # compile model
    model.compile(optimizer='adam', loss='mean_squared_error')
    ##output df
    output_df = pd.DataFrame(columns=['Trait', 'Task', 'RMSE', 'Correlation', 'Best_Prediction'])
    traits = list(["silk", "tassel", "oil", "protein", "starch"])
    tasks = list(["trait_mean", "trait_95_perc"])
    for i in traits:
        for j in tasks:
            model.fit(sim_x, sim_y, epochs=500, batch_size=1)
            predictions = model.predict(real_x)
            rmse = np.sqrt(mean_squared_error(real_y, predictions))
            correlation = np.corrcoef(real_y, predictions.flatten())[0, 1]
            best_predict = (np.argmax(array_A) == np.argmax(array_B))
            output_df = output_df.append({'Trait': trait, 'Task': task, 'RMSE': rmse,
                                          'Correlation': correlation, 'Best_Prediction': best_pred},
                                         ignore_index=True)