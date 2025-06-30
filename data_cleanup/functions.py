import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt



def show_PCA(feature_spaces: dict, regions: list):
    for region in regions:
        df = feature_spaces[region]
        antigen_labels = df["antigen_name"].tolist()

        # PCA (df darf keine numerische Spalten enthalten, außer Koordinaten))
        meta_spalten = ["pdb", "antigen_name", "antigen_species"]
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        X_2d = PCA(n_components=2).fit_transform(X)

        #### Farbzuordnung vorbereiten
        unique_antigens = sorted(df["antigen_name"].unique())
        color_map = plt.get_cmap('gist_rainbow_r', len(unique_antigens))
        color_dict = dict(zip(unique_antigens, [color_map(i) for i in range(len(unique_antigens))]))
        colors = [color_dict[a] for a in antigen_labels]


        color_map = plt.get_cmap('gist_rainbow_r', len(unique_antigens))
        color_dict = {antigen: color_map(i) for i, antigen in enumerate(unique_antigens)}
        colors = [color_dict[a] for a in antigen_labels]



        fig, axes = plt.subplots(2, 3, figsize=(12, 6))


        plt.figure(figsize=(8, 6))
        plt.scatter(X_2d[:, 0], X_2d[:, 1], c=pd.Categorical(feature_space["antigen_name"]).codes, cmap='tab10', s=10)
        plt.title(f"Embeddings der {region}")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.grid(True)

        plt.tight_layout()
        plt.show()



def show_feature_space(feature_space, CDR_region, antigen_labels):#zusätzliches Argument für color-coding: antigen
    
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
               for ant in unique_antigens]

    # Erzeuge Legende
    legend = legend_ax.legend(handles=handles, loc='center', ncol=1, fontsize='small')  # ncol=1 → alle untereinander

    legend_fig.tight_layout()
    plt.show()