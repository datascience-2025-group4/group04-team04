import os
os.environ["OMP_NUM_THREADS"] = "3"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines

from sklearn.decomposition import PCA

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

from collections import Counter
from scipy.stats import chi2_contingency

import math

regions = ["CDR_H1", "CDR_H2", "CDR_L1", "CDR_L2", "CDR_L3"]

##########
##########

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
        ax.set_title(f"PCA der {region}-Sequenzen")
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

##########
##########

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
        ax.set_xlabel('Anzahl der Cluster')
        ax.set_ylabel('Inertia (Summe der quadrierten Abstände)')
        ax.set_title(f'Elbow-Kurve für ({region})')
        ax.set_xticks(cluster_range) #x achse werden dann anzhal cluster gezeigt
        ax.grid(True) #gitternetzlinien im plot

    fig.tight_layout()
    fig.savefig("data/elbow_overview.png", dpi=150, bbox_inches='tight')
    plt.show()
        
##########
##########

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
        cluster_range = range(2, 20)  # Teste 2-10 Cluster

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
        ax.set_xlabel('Anzahl der Cluster')
        ax.set_ylabel('Silhouette Score')
        ax.set_title(f'Durchschnittliche Silhouette-Scores für ({region})')
        ax.set_xticks(cluster_range) #x achse werden dann anzhal cluster gezeigt
        ax.grid(True) #gitternetzlinien im plot

    fig.tight_layout()
    fig.savefig("data/silhouette_overview.png", dpi=150, bbox_inches='tight')
    plt.show()

##########
##########

def silhouette_check(feature_spaces: dict, region_ncluster: dict, ncols: int):
    '''
    Funktion, um bei einer gegebenen Clusterzahl die Silhouette-Scores der einzelnen Datenpunkte zu plotten
    Die Funktion nutzt noch für alle CDRs dieselbe Clusteranzahl auch, wenn unterschiedliche Clusterzahlen optimal wären!'''
    nrows = math.ceil(len(regions) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()
    
    for i, region in enumerate(regions):
        df = feature_spaces[region]

        n_clusters = region_ncluster[region]

        # PCA (ggf.)
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        #coords = PCA(n_components=2).fit_transform(X)
        coords = X

        # KMeans-Clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(coords)
    
        silhouette_avg = silhouette_score(coords, cluster_labels)
        sample_silhouette_values = silhouette_samples(coords, cluster_labels)

        # Silhouette-Scores aller Datenpunkte plotten
        ax = axes[i]
        y_lower = 10

        for j in range(n_clusters):
            silhouette_values = sample_silhouette_values[cluster_labels == j]
            silhouette_values.sort()
            cluster_size = silhouette_values.shape[0]
            y_upper = y_lower + cluster_size

            color = cm.nipy_spectral(float(j) / n_clusters)
            ax.fill_betweenx(
                np.arange(y_lower, y_upper),
                0,
                silhouette_values,
                facecolor=color,
                edgecolor=color,
                alpha=0.7
            )

            ax.text(-0.05, y_lower + 0.5 * cluster_size, str(j))
            y_lower = y_upper + 10

        ax.axvline(x=silhouette_avg, color="red", linestyle="--")
        ax.set_title(f"Silhouette Scores für optimale\nCluster-Anzahl für {region}")
        ax.set_xlabel("The silhouette coefficient values")
        ax.set_ylabel("Cluster label")
        ax.set_yticks([])
        ax.set_xlim([-0.1, 1.0])
        ax.set_ylim([0, len(coords) + (n_clusters + 1) * 10])
    
    fig.tight_layout()
    fig.savefig("data/check_silhouette.png", dpi=150, bbox_inches="tight")
    plt.show()

##########
##########

def show_Kmeans(feature_spaces: dict, region_ncluster: dict, n_cols: int):
    '''
    Clustering auf basis der PCA
    kmeans mit 2 => laut elbow method am sinnvollsten
    '''
    nrows = math.ceil(len(regions) / n_cols)
    fig, axes = plt.subplots(nrows, n_cols, figsize=(5 * n_cols, 4.5 * nrows))
    axes = axes.flatten()
    
    for i, region in enumerate(regions):
        df = feature_spaces[region]
        n_clusters = region_ncluster[region]
        antigen_labels = df["antigen_name"].tolist()


        # PCA (ggf.)
        koordinaten_spalten = df.select_dtypes(include='number').columns
        X = df[koordinaten_spalten].values
        #coords = PCA(n_components=2).fit_transform(X)
        coords = X

        # KMeans-Clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(coords)

        ax = axes[i]
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=cluster_labels, cmap='tab10', s=20, alpha=0.7)
        ax.set_title(f"PCA + KMeans Clustering\n({region}, k={n_clusters})")
        ax.set_xlabel("PCA 1")
        ax.set_ylabel("PCA 2")
        ax.grid(True)

        # Farbskala für die Cluster-Farben
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label("Cluster Label")

        # Antigen-Häufigkeiten in den Clustern
        # Erzeuge Dictionary: cluster_label --> Liste von Antigen-Names
        cluster_antigen_mapping = {}
        for cluster_label, antigen_name in zip(cluster_labels, antigen_labels):
            # Wenn cluster-Label noch nicht im Dictionary: leeres Array initialisieren
            if cluster_label not in cluster_antigen_mapping:
                cluster_antigen_mapping[cluster_label] = []
            # Füge den Antigen-Namen zur Liste für dieses Cluster hinzu
            cluster_antigen_mapping[cluster_label].append(antigen_name)
        # Jetzt haben wir für jedes Cluster eine Liste mit den vorkommenden Antigen-Namen

        print(f"Antigen-Häufigkeiten für {region} (k={n_clusters})")
        for cluster_label in sorted(cluster_antigen_mapping):
            print(f"\nCluster {cluster_label}:")
            for antigen, count in Counter(cluster_antigen_mapping[cluster_label]).most_common():
                print(f"    {antigen}: {count}x")
        
            # Häufigkeit der Antigen-Namen im Cluster
            counter = Counter(cluster_antigen_mapping[cluster_label])
            # Antigen-Namen und Anzahl ausgeben, sortiert nach Häufigkeit
            for antigen, count in counter.most_common():
                print(f"{antigen}: {count}")

    fig.tight_layout()
    fig.savefig("data/kmeans_clustering.png", dpi=150, bbox_inches="tight")
    plt.show()

##########
##########

def chi_sq():
    '''
    ...
    '''
    #Proportionstest => testen, ob die Antigen-Namen systematisch unterschiedlich auf die Cluster verteilt sind

    #nullhypothese (H0): Verteilung der Antigen-Namen auf die Cluster ist zufällig
    #H1: Verteilung der Antigen-Namen auf die Cluster ist nicht zufällig

    for region in regions:

        # Kontingenztabelle (Antigen x Cluster) erzeugen

        # automatische Zählung und erstellen von dataframe über crosstab
        contingency_table = pd.crosstab(pd.Series(antigen_labels, name='Antigen'), pd.Series(cluster_labels, name='Cluster'))

        #kontingenztabelle printen
        print(f"\n Kontingenztabelle für {region}")
        print(contingency_table)

        #chi-quadrat-test
        chi2, p, dof, expected = chi2_contingency(contingency_table)
        #chi2=maß für wie stark die beobachteten häufigkeiten von den erwarteten abweichen (je größer, desto größer abweichung)
        #p-Wert=Wahrscheinlichkeit, bei zufälliger Verteilung eine Abweichung mindestens so groß wie chi2 zu beobachten
        #dof=degrees of freedom (#zeilen-1 * #spalten-1)
        #expected=tabelle mit erwarteten häufigkeiten unter nullhypothese

        #print ergebnisse
        print(f"\n Chi-Quadrat-Test Ergebnis für {region}")
        print(f"Chi2-Wert: {chi2:.4f}") #4f für 4 nachkommSTELLEN
        print(f"p-Wert:    {p:.4e}")  # 4e für exponentielle darstellung mit 4 nachkommastellen
        print(f"Freiheitsgrade: {dof}")

        #print die expected werte zum vergleich
        expected_df = pd.DataFrame(expected, index=contingency_table.index, columns=contingency_table.columns)
        print(f"\nErwartete Häufigkeiten (unter Nullhypothese):")
        print(expected_df.round(2)) #auf 2 nachkommastellen runden

        #die tatsache dass die antigen_name häufigkeiten sehr unterschiedlich sind wird hier automatisch berücksichtigt: 
        #expected-Array wird aus den Randhäufigkeiten berechnet => randhäufigkeiten für antigen =  Summe pro Zeile => wie oft kommt jedes Antigen insgesamt im ganzen Datensatz vor?
        #randhäufigkeiten für cluster=Summe pro Spalte => wie groß ist jedes Cluster insgesamt?

        #wenn p < 0.05 => signifikant also die verteilung unterschiedet sich signifikant von zufall => N0 hypothese muss verworfen werden denn die unterschiede in der verteilung sind nicht durch zufall zu erklären
