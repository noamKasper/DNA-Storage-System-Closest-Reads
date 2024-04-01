import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from brokenaxes import brokenaxes
from n_closest_post import POST_PROCESSING_FULL_RESULTS_PATH

N_CLOSEST_FIG_FOLDER = "/home/noam/alphaProject/figs/n_closest"

def bad_fig(df, path):
    df.reset_index(inplace=True)
    df = pd.melt(df, id_vars="dataset", var_name='classification type', value_name='count')
    df = df.iloc[4:, :]

    print(df)
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=df, x="dataset", y="count", hue="classification type")
    ax.set_yscale('log')  # Set y-axis to logarithmic scale
    plt.xlabel('Database')
    plt.ylabel('Count (log scale)')
    plt.title('Classification Count per Database N=5')
    plt.savefig(path)


def improved_stack_plot(df: pd.DataFrame, path, ascending=True, text_size=12, text_size_amp=2):
    databases = df['dataset'].unique()  # Get unique database names
    df.drop('dataset', axis=1, inplace=True)

    # Normalize data to ensure each bar reaches 100%
    df_normalized = df.div(df.sum(axis=1), axis=0) * 100

    # Custom color mapping
    colors = {'LT': '#779ecc', 'ST': '#9fc0de', 'LF': '#ff985a', 'SF': '#f2c894'}

    plt.figure(figsize=(12, 8))
    bax = brokenaxes(ylims=((0, 10.5), (97.5, 100)), hspace=0.03)
    bottom = pd.Series([0] * len(df_normalized))

    # Sort columns by sum and reverse if needed
    sorted_columns = df_normalized.sum().sort_values(ascending=ascending).index

    for col in sorted_columns:
        bax.bar(databases, df_normalized[col], bottom=bottom,
                label=col, color=colors[col], width=0.5, edgecolor='black', linewidth=1)  # Rounded edges and thicker outline
        bottom += df_normalized[col]

    # Add labels and title
    bax.set_xlabel('\nDatabase', fontsize=text_size + text_size_amp)  # Increase font size
    bax.set_ylabel('Percentage', fontsize=text_size + text_size_amp)  # Increase font size
    bax.set_title('Classification Count Percentage per Database N=5\n', fontsize=text_size + text_size_amp * 2)  # Set title font size
    bax.tick_params(axis='both', which='major', labelsize=text_size)  # Increase tick label font size
    bax.legend(fontsize=text_size, loc='best', reverse=True).set_zorder(1.5)  # Set legend font size
    bax.locator_params(axis='x', nbins=len(databases))


    plt.savefig(path)


if __name__ == "__main__":
    df = pd.read_csv(POST_PROCESSING_FULL_RESULTS_PATH)
    improved_stack_plot(df, os.path.join(N_CLOSEST_FIG_FOLDER, "5.png"), text_size=13)

