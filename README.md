# group04-team04
See the project description here: https://github.com/maiwen-ch/2025_Data_Analysis_Topic_04_Antibody_Antigen_Interactions

<details>
<summary><strong>Contributors</strong></summary>

Marcel Buß, Pauline Kipp, Christoph Plenge, Anna van der Ham
</details>

<details>
<summary><strong>Supervisor</strong></summary>

Prof. Dominik Niopek ([dominik.niopek@uni-heidelberg.de])  
Benedict Wolf ([benedict.wolf@uni-heidelberg.de])
Jan Mathony ([jan.mathony@uni-heidelberg.de]) 

Tutor: Enno Schaefer ([enno.schaefer@uni-heidelberg.de])   
</details>

<details>
<summary><strong>Project Goal</strong></summary>
The goal of our project was to investigate whether the established Canonical Forms of antibody CDRs can be reconstructed using different computational approaches. We explored a naive method based on residue percentages, as well as more advanced techniques such as ESM embeddings and MMseqs2 clustering.

Furthermore, we analyzed how different antigens are distributed across these canonical forms, aiming to identify potential structural or functional preferences of antigens for specific CDR conformations.
</details>

<details>
<summary><strong>Structure of the repository</strong></summary>
Before running the code, please make sure to download the folder "data" from the link provided at the end of this document (Additional files and folders). Additionally, data must be placed directly inside the group04_team04 directory. If the folder is renamed or placed elsewhere, the code will not function correctly.






In our repository, the final notebook, that generates all important plots is found in the Documentation folder and is called [P53_DMS_Documentation.ipynb](Documentation/P53_DMS_Documentation.ipynb).
The [report (as a pdf)](Documentation/report_DMS_topic02_team02.pdf) can be found in the same folder. Within our Documentation, you will find **five sub-topics**. The first one looks at the Comparability of p53 Datasets.
The code generating the relevant plots can be found [within the Documentation folder](Documentation/backgrounddata.py) and 
for visualization purposes, that might also be used later on, code from the [Visualization folder](visualization) was used.
The code for plots on the other four topics can be found here: 
- [Data cleanup](data_cleanup) 
- [Data exploration](data_exploration)
- [Domain comparison](domain_comparison) 
- [Calculating severity scores](severity_score)

In each of these folders, the relevant functions are defined in the .py file. In most folders, there is an additional **exploratory**
folder containing all the experimental notebooks. Jupyter notebooks, that are within the sub-topic folders but not in the exploratory 
folders are mentioned in the [P53_DMS_Documentation.ipynb file](Documentation/P53_DMS_Documentation.ipynb) and contain further 
information, outlook or important additional information. 

</details>




<details>
<summary><strong>Covering the mandatory aspects of the project</strong></summary>

Our project was supposed to contain the following elements. Here, we list which sub-topic covers which mandatory aspect: 
- **descriptive statistics** analyzing antigen distribution in the original dataset (data exploration) to evalute vadility of Chi2 Test
- **graphical representations**: (data_exploration, naiver_ansatz, esmc_ansatz, proportionstest, v_measure)
- **dimension reduction** PCA (esmc_ansatz), hierarchical clustering (naiver_ansatz, esmc_ansatz), MMseqs2 Clustering (MMseqs2)
- **statistical tests** Chi2 Test of Independencet (Proportionstest)
- **linear regression** was not implemented (as discussed with our tutor Enno Schäfer)

</details>

<details>
<summary><strong>Additional files and folders</strong></summary>
-environment file
-data folder: link

</details>

<details>
<summary><strong>Download the datasets worked on</strong></summary>


</details>