import numpy as np
import pandas as pd
import os
import torch
# torch ist ein Paket aus der library PyTorch zur Arbeit mit neuronalen Netzen
from torch import nn
# nn ist ein Paket aus Torch, mit dem man neuronale Netze bauen kann
from esm.models.esmc import ESMC
from esm.utils.constants.models import ESMC_600M
from esm import pretrained
from esm.tokenization import EsmSequenceTokenizer
# EsmSequenceTokenizer ist eine Objektklasse, Funktionen wir .encode und .decode hat,
# um zwischen Sequenz und Aminosäure-Tokens zu unterscheiden
from esm.utils import encoding
# encoding ist ein Modul (.py-Datei), das Funktionen wie tokenize_sequence(seq=str, tokenizer=EsmSequenceTokenizer, ...) enthält
from torch.cuda.amp import autocast
# wichtig für Ausführung des Codes auf der GPU

esm = ESMC.from_pretrained("esmc_600m").to('cuda')
tokenizer = EsmSequenceTokenizer()

df_list = ['/gpfs/lsdf02/sd24c002/ProDomino_Extensions/tobi/data/TED/single_no_prodomino_caths_atall.csv',
          '/gpfs/lsdf02/sd24c002/ProDomino_Extensions/tobi/data/TED/single_very_strict.csv',
          '/gpfs/lsdf02/sd24c002/ProDomino_Extensions/tobi/data/TED/single_minus_prodomino.csv',
          '/gpfs/lsdf02/sd24c002/ProDomino_Extensions/tobi/data/TED/single_ready.csv']

df = pd.concat((pd.read_csv(f) for f in df_list), ignore_index=True)
df_no_dupes = (
    df
    .drop_duplicates(subset=['prot_id'], keep='first')
    .reset_index(drop=True)
)

to_process = [
    (n, s)
    for n, s in zip(df_no_dupes.prot_id, df_no_dupes.seq_no_insert)
    if not os.path.exists(f'…/Target_{n}.npy')
]

print(len(to_process))

with torch.no_grad():
# torch ist eine library, um mit neuronalen Netzen zu arbeiten
# no_grad heißt: keine Anpassung der weights und bias um loss-function zu minimieren, sondern neuronales Netz nur "vorwärts" laufen lassen,
# weil wir das Modell nutzen und nicht trainieren wollen               
    for name, seq in to_process:
        tokenized = encoding.tokenize_sequence(seq, tokenizer, add_special_tokens=True).to('cuda')
        # encoding ist ein Modul aus dem Paket esm.utils
        pred = esm.forward(tokenized.unsqueeze(0))
        embeddings = pred.embeddings.to('cpu').squeeze()
        as_arr = embeddings.float().detach().numpy()
        as_arr = as_arr[1:-1,:]
        np.save(f'/gpfs/lsdf02/sd24c002/ProDomino_Extensions/data/TED_ESMC600/Target_{name}.npy', as_arr)
