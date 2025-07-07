import os
os.environ["OMP_NUM_THREADS"] = "3"
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.decomposition import PCA


def show_PCA(feature_spaces: dict, seq_regions : list):
    '''
    Input:
    - feature_spaces: Dictionary, enthält für jede CDR-Region einen Feature-Space
    --> Feature-Space hat folgende Struktur: ["pdb" + "antigen_name" + "antigen_species" + Koordinaten]
    - seq_regions = ["SEQ_H1", "SEQ_H2", "SEQ_L1", "SEQ_L2", "SEQ_L3"]
    Output:
    - Visualisierung der 2D-PCAs aller CDR-Regionen
    '''
    ncols = 3
    nrows = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    for i, seq_region in enumerate(seq_regions):
        df = feature_spaces[seq_region]
        
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
        ax.set_title(f"PCA der {seq_region}-Sequenzen")
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


def show_single_PCA(feature_space, region, method):
    '''
    Input:
    - feature_spaces: Dictionary, enthält für jede CDR-Region einen Feature-Space
    --> Feature-Space hat folgende Struktur: ["pdb" + "antigen_name" + "antigen_species" + Koordinaten]
    - region: "SEQ_H1", or "SEQ_H2", or "SEQ_L1", or "SEQ_L2", or "SEQ_L3"
    - method: "ESMC", or "frequency"
    Output:
    - Visualisierung der 2D-PCAs aller CDR-Regionen
    '''
                   
    # PCA
    koordinaten_spalten = feature_space.select_dtypes(include='number').columns
    X = feature_space[koordinaten_spalten].values
    coords = PCA(n_components=2).fit_transform(X)

    # Farbzuordnung vorbereiten
    antigen_names = sorted(feature_space["antigen_name"].unique())
    color_map = plt.get_cmap('gist_rainbow_r', len(antigen_names))
    color_dict = {}
    for j, antigen in enumerate(antigen_names):
        color = color_map(j)
        color_dict[antigen] = color
    colors = []
    for antigen in feature_space["antigen_name"]:
        colors.append(color_dict[antigen])

    # PCA plotten mit Color-Coding
    plt.scatter(coords[:, 0], coords[:, 1], c=colors, s=20, alpha=0.7)
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.title(f"PCA der {region}-Sequenzen")
    plt.grid(True)

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
    plt.legend(
        handles = handles,
        fontsize = 'xx-small',
        loc = 'best',
        framealpha = 0.3,
        markerscale = 0.7,
        labelspacing = 0.2
    )
    plt.tight_layout()
    plt.savefig(f"data/pca_{region}_{method}.png", dpi=150, bbox_inches='tight')
    plt.show()