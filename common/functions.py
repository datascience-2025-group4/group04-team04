def Notebook_vorbereiten():
    """ Diese Funktion bereitet das Notebook vor, indem es Standardbibliotheken importiert und die Rohdaten lädt."""
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    ab_ag = pd.read_csv('download/ab_ag.tsv', sep='\t')
    columns = pd.read_csv('download/columns.tsv', sep='\t')
    uniprot_data = pd.read_csv('download/uniprot_data.tsv', sep='\t')