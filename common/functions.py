def Notebook_vorbereiten():
    """ Diese Funktion bereitet das Notebook vor, indem es Standardbibliotheken importiert und die Rohdaten lädt."""
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    ab_ag_raw = pd.read_csv('data/ab_ag.tsv', sep='\t')
    columns = pd.read_csv('data/columns.tsv', sep='\t')
    uniprot_data = pd.read_csv('data/uniprot_data.tsv', sep='\t')
    return ab_ag_raw, columns, uniprot_data