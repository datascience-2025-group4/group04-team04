import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math

regions = ["CDR_H1", "CDR_H2", "CDR_H3", "CDR_L1", "CDR_L2", "CDR_L3"]


def show_PCA(feature_spaces: dict, ncols: int):
    '''
    Input:
    - feature_spaces: Dictionary, enthält für jede CDR-Region einen Feature-Space
    - Feature-Space hat folgende Struktur: ["pdb" + "antigen_name" + "antigen_species" + Koordinaten]
    Output:
    - Visualisierung der 2D-PCAs aller CDR-Regionen
    '''
    fig, axes = plt.subplots(1, 6, figsize=(8, 4 * len(regions)))
    nrows = math.ceil(len(regions) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    for i, region in enumerate(regions):
        df = feature_spaces[region]
        
        # PCA
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        coords = PCA(n_components=2).fit_transform(X)

        # Farbzuordnung vorbereiten
        antigen_names = sorted(df["antigen_name"].unique())
        color_map = plt.get_cmap('gist_rainbow_r', len(antigen_names))
        color_dict = {}
        for j, antigen in enumerate(antigen_names):
            color = color_map(j)
            color_dict[antigen] = color
        colors = []
        for antigen in df["antigen_name"]:
            colors.append(color_dict[antigen])

        # PCA plotten mit Color-Coding
        ax = axes[i]
        ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=20, alpha=0.7)
        ax.set_xlabel("PCA 1")
        ax.set_ylabel("PCA 2")
        ax.set_title(f"PCA of {region} Sequences")
        ax.grid(True)
        

        # Legende zum Color-Coding
        handles = []
        for antigen in antigen_names:
            color = color_dict[antigen]
            handle = mlines.Line2D(
                [], [],
                color = color,
                marker = "o",
                linestyle = "None",
                markersize = 6,
                label = antigen 
            )
            handles.append(handle)
        ax.legend(
            handles = handles,
            fontsize = 'xx-small',
            loc = 'best',
            framealpha = 0.3,
            markerscale = 0.7,
            labelspacing = 0.2
        )
    fig.tight_layout()
    fig.savefig("data/pca_overview.png", dpi=150, bbox_inches='tight')
    plt.show()