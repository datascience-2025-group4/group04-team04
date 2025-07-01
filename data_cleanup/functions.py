import os
os.environ["OMP_NUM_THREADS"] = "3"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines

from sklearn.decomposition import PCA

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

from collections import Counter
from scipy.stats import chi2_contingency

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


def elbow(feature_spaces: dict, ncols: int):
    '''
    Input:
    - feature_spaces: Dictionary, enthält für jede CDR-Region einen Feature-Space
    - Feature-Space hat folgende Struktur: ["pdb" + "antigen_name" + "antigen_species" + Koordinaten]
    Output:
    - Visualisierung der Elbow-Kurven, um optimale Cluster-Zahl für K-Means zu bestimmen
    '''
    nrows = math.ceil(len(regions) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    for i, region in enumerate(regions):
        df = feature_spaces[region]
        
        # PCA (ggf.)
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        #coords = PCA(n_components=2).fit_transform(X)
        coords = X

        # Range of cluster numbers to try
        cluster_range = range(1, 11)  # Teste 1-10 Cluster

        # Hier Inertia-Werte speichern, [] für liste
        #inertia => Summe der quadrierten Abstände der Punkte zum jeweiligen Cluster-Zentrum => je kleiner, desto besser
        inertia_values = []

        # Berechne KMeans für jede Cluster-Anzahl
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            kmeans.fit(coords)  # fit passt daten an modell an => kmeans bekommt daten (coords der antikörper punkte) und initialisiert Cluster Zentren und macht dann iterative optimierung 
            # also weist punkte dem nächsten cluster zentrum zu, bei den gebildeten cluster wenn neue zentren berechnet und die punkte werden neu zugeordnet
            inertia_values.append(kmeans.inertia_)

        # Elbow-Kurven plotten
        ax = axes[i]
        ax.plot(cluster_range, inertia_values, marker='o')
        ax.set_xlabel('Anzahl der Cluster (n_clusters)')
        ax.set_ylabel('Inertia (Summe der quadrierten Abstände)')
        ax.set_title(f'Elbow-Methode zur Bestimmung optimaler Clusterzahl ({region})')
        ax.set_xticks(cluster_range) #x achse werden dann anzhal cluster gezeigt
        ax.grid(True) #gitternetzlinien im plot

    fig.tight_layout()
    fig.savefig("data/elbow_overview.png", dpi=150, bbox_inches='tight')
    plt.show()
        


def silhouette(feature_spaces: dict, ncols: int):
    '''
    ...
    '''
    nrows = math.ceil(len(regions) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    for i, region in enumerate(regions):
        df = feature_spaces[region]
        
        # PCA (ggf.)
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        #coords = PCA(n_components=2).fit_transform(X)
        coords = X

        # Range of cluster numbers to try
        cluster_range = range(1, 11)  # Teste 1-10 Cluster

        silhouette_values = []

        # Berechne Silhouette Score für jede Clusteranzahl
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            labels = kmeans.fit_predict(coords)
            score = silhouette_score(coords, labels)
            silhouette_values.append(score)           

        # Silhouett-Werte plotten
        ax = axes[i]
        ax.plot(cluster_range, silhouette_values, marker='o')
        ax.set_xlabel('Anzahl der Cluster (n_clusters)')
        ax.set_ylabel('Silhouette Score')
        ax.set_title(f'silhouette-Methode zur Bestimmung optimaler Clusterzahl ({region})')
        ax.set_xticks(cluster_range) #x achse werden dann anzhal cluster gezeigt
        ax.grid(True) #gitternetzlinien im plot

    fig.tight_layout()
    fig.savefig("data/silhouette_overview.png", dpi=150, bbox_inches='tight')
    plt.show()


def silhouette_check(feature_spaces: dict, ncluster: int, ncols: int):
    '''
    Funktion, um bei einer gegebenen Clusterzahl die Silhouette-Scores der einzelnen Datenpunkte zu plotten'''
    nrows = math.ceil(len(regions) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()
    
    for i, region in enumerate(regions):
        df = feature_spaces[region]

        # PCA (ggf.)
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        #coords = PCA(n_components=2).fit_transform(X)
        coords = X

        kmeans = KMeans(n_clusters=ncluster, random_state=42)
        cluster_labels = kmeans.fit_predict(coords)
    
        silhouette_avg = silhouette_score(coords, cluster_labels)
        sample_silhouette_values = silhouette_samples(coords, cluster_labels)

        # Plot
        fig, ax1 = plt.subplots(1, 1)
        fig.set_size_inches(8, 6)

        y_lower = 10
        for i in range(n_clusters):
            ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.nipy_spectral(float(i) / n_clusters)
            ax1.fill_betweenx(
                np.arange(y_lower, y_upper),
                0,
                ith_cluster_silhouette_values,
                facecolor=color,
                edgecolor=color,
                alpha=0.7
            )

            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
            y_lower = y_upper + 10

        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
        ax1.set_title(f"Silhouette analysis for KMeans clustering: {region}")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
        ax1.set_yticks([])
        ax1.set_xlim([-0.1, 1.0])
        ax1.set_ylim([0, len(coords) + (n_clusters + 1) * 10])
        plt.tight_layout()
        plt.show()
def show_Kmeans():
    '''
    ...
    '''

def chi_sq():
    '''
    ...
    '''