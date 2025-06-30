import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

regions = ["CDR_H1", "CDR_H2", "CDR_L1", "CDR_L2", "CDR_L3"]


def show_PCA(feature_spaces: dict):
    '''
    Input:
    - feature_spaces: Dictionary, enthält für jede CDR-Region einen Feature-Space
    - Feature-Space hat folgende Struktur: ["pdb" + "antigen_name" + "antigen_species" + Koordinaten]
    Output:
    - Visualisierung der 2D-PCAs aller CDR-Regionen
    '''
    fig, axes = plt.subplots(5, 0, figsize=(15, 8))
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
        ax.legend(handles=handles, fontsize='x-small', loc='best')
    fig.tight_layout()
    plt.show()



def show_feature_space(feature_space, CDR_region, antigen_names):#zusätzliches Argument für color-coding: antigen
    
    # Plotten der PCA-Koordinaten in 2D
    fig, (ax_scatter, ax_legend) = plt.subplots(2,1,figsize=(10,10), gridspec_kw={'height_ratios': [4,1]})
    ax_scatter.scatter(coords[:, 0], coords[:, 1], c=colors, s=20, alpha=0.7)#zusätzliche Argumente für color-coding: c=color_codes, cmap='tab20', 
    ax_scatter.set_xlabel("PCA 1")
    ax_scatter.set_ylabel("PCA 2")
    ax_scatter.set_title(f"PCA of {CDR_region} Sequences")


    # 2. Legende als eigene Figure
    legend_fig, legend_ax = plt.subplots(figsize=(12, 12))  # Größe anpassen je nach Anzahl Antigene
    legend_ax.axis('off')  # Keine Achsen

    # Erstelle Handles
    handles = [mlines.Line2D([], [], color=color_dict[ant], marker='o', linestyle='None',
                             markersize=8, label=ant)
               for ant in antigen_names]

    # Erzeuge Legende
    legend = legend_ax.legend(handles=handles, loc='center', ncol=1, fontsize='small')  # ncol=1 → alle untereinander

    legend_fig.tight_layout()
    plt.show()