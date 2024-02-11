import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from n_closest_post import POST_PROCESSING_FULL_RESULTS_PATH

N_CLOSEST_FIG_FOLDER = "/home/noam/alphaProject/figs/n_closest/5.png"

if __name__ == "__main__":
    df = pd.read_csv(POST_PROCESSING_FULL_RESULTS_PATH)
    df.reset_index(inplace=True)
    df = pd.melt(df, id_vars="dataset", var_name='classification type', value_name='count')
    df = df.iloc[4:, :]

    print(df)
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=df, x="dataset", y="count", hue="classification type")
    ax.set_yscale('log')  # Set y-axis to logarithmic scale
    plt.title('Classification Count per Database')
    plt.savefig(N_CLOSEST_FIG_FOLDER)

