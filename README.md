# group04-team04
See the project description here: https://github.com/maiwen-ch/2025_Data_Analysis_Topic_04_Antibody_Antigen_Interactions

<details>
<summary><strong>Contributors</strong></summary>

Marcel Buß, Pauline Kipp, Christoph Plenge, Anna van der Ham
</details>

<details>
<summary><strong>Supervisor</strong></summary>

Prof. Dr. Dominik Niopek ([dominik.niopek@uni-heidelberg.de]) 
Dr. Jan Mathony ([jan.mathony@uni-heidelberg.de])  
Benedict Wolf ([benedict.wolf@uni-heidelberg.de])


Tutor: Enno Schaefer ([enno.schaefer@uni-heidelberg.de])   
</details>

<details>
<summary><strong>Project Goal</strong></summary>
The goal of our project was to investigate whether the established Canonical Forms of antibody CDRs can be reconstructed using different computational approaches. We explored a naive method based on residue percentages, as well as more advanced techniques such as ESM Cambrian embeddings and MMseqs2 clustering.

Furthermore, we analyzed how different antigens are distributed across these canonical forms, aiming to identify potential structural or functional preferences of antigens for specific CDR conformations.
</details>

<details>
<summary><strong>Structure of the repository</strong></summary>

The repository is organized into the following key directories.
The data_cleanup folder[data_cleanup] filters all relevant rows from the raw data and then uses SCALOP to extract both canonical forms and the CDR regions. The data_exploration folder [data_exploration] builds visualizations of antigen distribution in the filtered datset. There are three different folder for our three approaches: naive method [naiver_ansatz], ESMC [esmc_ansatz], MMseqs2 [MMseqs2]. The v_measure directory [v_measure] compares our own clustering results against reference clusters to quantify their agreement, for our three approaches. In the folder chi2 [chi2], chi-square tests are applied to analyze antigen distributions within different clusters. 

Finally, the data folder [data] houses all project-generated files. Before running the code, please make sure to download the folder data from the link provided at the end of this document (Additional files and folders). Additionally, data must be placed directly inside the group04_team04 directory. If the folder is renamed or placed elsewhere, the code will not function correctly. 

The folder data already contains all files which are created if the code is run again

In our repository the final plots, v measure scores or p-values are generated within the corresponding folders.

</details>

<details>
<summary><strong>Covering the mandatory aspects of the project</strong></summary>

Our project was supposed to contain the following elements. Here, we list which sub-topic covers which mandatory aspect: 
- **descriptive statistics** analyzing antigen distribution in the original dataset (data exploration) to evalute vadility of Chi2 Test
- **graphical representations**: (data_exploration, naiver_ansatz, esmc_ansatz, proportionstest, v_measure)
- **dimension reduction** PCA (esmc_ansatz), hierarchical clustering (naiver_ansatz, esmc_ansatz), MMseqs2 Clustering (MMseqs2)
- **statistical tests** Chi2 Test of Independencet (Proportionstest)
- **linear regression** was not implemented (as discussed with our tutor Enno Schäfer)

<details>
<summary><strong>requirements</strong></summary>

for embedding.ipynb: GPU 
for MMseqs2 and assign function from SCALOP: LinuxOS


<details>
<summary><strong>Additional files and folders</strong></summary>

-the environment file [environment.yml] contains all packages to run the code
-data folder:- https://heibox.uni-heidelberg.de/d/308698a78f2043c292bf/

SAbDab download script:
- used in [data_cleanup/data_cleanup.ipynb] to download VH and VL sequences from the SAbDab
- https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/downloads/sabdab_downloader.py/
- already included in [data] folder, if our repository is executed as advised

SCALOP:
- used in [data_cleanup/canonical_forms_extrahieren.ipynb] to extract CDR sequences and canonical forms
- https://github.com/oxpig/SCALOP.git 
- relevant packages and functions already included in [environment.yml]

MMseqs2_
- used in [MMseqs2/MMseqs2.ipynb] for clustering of CDR sequences
- https://github.com/soedinglab/mmseqs2/wiki#getting-started
- please follow the installation steps in the README file of the MMseqs2 repository and execute in LinuxOS or WSL terminal
</details>