## Importation des fichiers edf contenant les signaux de PSG

## Remarque : Le format .edf permet de rassembler des signaux n'ayant pas la même fréquence d'échantillonnage.

## Lecture du fichier de PSG grâce au package .mne
## Biological time series : MNE library an open-source Python module for processing, analysis, and visualization of functional neuroimaging data (https://mne.tools/dev/install/index.html)

import os
import txt
import pandas as pd
import numpy as np
import re
import striprtf
from striprtf.striprtf import rtf_to_text
import datetime as dt
import glob
import mne
import matplotlib.pyplot as plt

# Fonction qui
    ## vérifie qu'il existe un fichier .edf dans le dossier indiqué;
    ## renvoie un message d'erreur en cas d'absence de fichier sinon le fichier data edf chargé
def verifier_isPSG(path, patient):
    try :
        file = glob.glob(path + patient + "/*.edf")[0]
    except : # Exception : il n'y a pas de fichier edf dans le dossier.
        raise Exception("Pas de fichier EDF")
    else :
        data = mne.io.read_raw_edf(file, preload=True) # Chargement des données.
        return data

## Récupération des signaux d'intérêt

## SATURATION EN OXYGÈNE
## Filtrage du signal de saturation
def filtrer_signalSat(sat_df):
    # Filtrage des Valeurs non physiologiques
    sat_df = sat_df.iloc[np.where((sat_df.SAT > 60) & (sat_df.SAT < 100))[0]]
    ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(sat_df)-1, num = len(sat_df), dtype = int)
    sat_df = sat_df.assign(index = index)
    sat_df.set_index(index, drop = True, inplace = True)
    sat_df.drop(columns = "index", inplace = True)
    # Filtrage des variations non physiologiques
    idx_drop5 = np.where(abs(np.diff(sat_df.SAT))>=5)[0] ## Repérage des indices avec drops > 5 et filtrage des drops
    mean = sat_df.SAT.mean() ## moyenne
    var = 5
    #np.std(sat_df.SAT) ## Ecart-type sur l'ensemble des données
    masque = []
    for i in idx_drop5 :

        if (i>= len(sat_df)):
            break

        masque.extend([i]) ## filtrage des drops
        k = 1
    ## filtrage de tous les points adjacents au drop hors de l'intervalle [mean - std : mean + std]
        while ((sat_df.loc[i + k, "SAT"] < mean - var) | (mean + var < sat_df.loc[i + k, "SAT"])): ## filtrage en aval du drop
            masque.extend([i+k])
            k += 1
            if (i+k >= len(sat_df)) :
                break
        k = 1
        while ((sat_df.loc[i - k, "SAT"] < mean - var) | (mean + var < sat_df.loc[i - k, "SAT"])): ## filtrage en amont du drop
            masque.extend([i-k])
            k += 1
            if (i-k <= 0) :
                break

    sat_df = sat_df.drop(masque)
    ## Mise en place d'une moyenne mobile pour terminer le filtrage des variations non physiologiques
        ## A REPRENDRE (entraine des problèmes dans le repérage des désaturations car enlève les 4 premiers points du signal)
    ## Idée : mettre en place un filtrage des outliers par un Modified Thompson Tau Test
        ## sur l'ensemble du signal ou par fenêtre de 5 secondes
    ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(sat_df)-1, num = len(sat_df), dtype = int)
    sat_df = sat_df.assign(index = index)
    sat_df.set_index(index, drop = True, inplace = True)
    sat_df.drop(columns = "index", inplace = True)

    return sat_df

# Récupération du signal correspondant à la saturation sous forme d'un DataFrame de 2 colonnes : Heure et SAT.
## Remarque : En cas d'ajout d'un nouveau logiciel, il faut ajouter le nom du signal correspondant à la saturation en oxygène dans le tableau "saturation"
def get_signalSat(data, annot_df):
    raw_data = data.get_data()
    channels = data.ch_names
    edf_df = pd.DataFrame(raw_data, index = channels)
    ## Récupération du signal de saturation.
    saturation = ["Sp02", "SAT"] ## Deltamed = SAT - Remlogic = "Sp02"
    for ch in channels :
        if (ch in saturation) :
            sat_df = pd.DataFrame({"SAT" : edf_df.loc[ch]})
            sfreq = data._raw_extras[0]["n_samps"][channels.index(ch)]
    ## Mise en place d'une colonne du temps en fonction du fichier d'annotation.
    sat_df = pd.DataFrame(sat_df.index, columns = ["Heure"]).join(sat_df, how = 'inner')
    sat_df["Heure"] = sat_df["Heure"]/sfreq
    heureDebut = annot_df.loc[0, "Heure"] # Récupération de l'heure, en seconde de départ
    sat_df["Heure"] += heureDebut # Modification de la colonne de temps du fichier .EDF

    sat_df = filtrer_signalSat(sat_df) ## Filtrage

    ## Optimisation de la mémoire
    sat_df.Heure = np.float32(sat_df.Heure)
    sat_df.SAT = np.float32(sat_df.SAT)
    sat_df.index = pd.RangeIndex(0,len(sat_df),1, dtype=np.int32)

    return (sat_df, sfreq)

## PPG
# Récupération du signal correspondant à la PPG sous forme d'un DataFrame normé.
def get_signalPPG(data, annot_df):
    raw_data = data.get_data()
    channels = data.ch_names
    edf_df = pd.DataFrame(raw_data, index = channels)
## Récupération du signal de PPG.
    PPG = ["Pouls", "Plethysmogram"] ## Deltamed = Pouls - Remlogic = "Plethysmogram"
    for n in range(0, len(channels)):
        if (channels[n] in PPG):
            PPG_df = pd.DataFrame({"PPG" : edf_df.loc[channels[n]]})
            sfreq = data._raw_extras[0]["n_samps"][n]
## Mise en place d'une colonne du temps en fonction du fichier d'annotation.
    PPG_df = pd.DataFrame(PPG_df.index, columns = ["Heure"]).join(PPG_df, how = 'inner')
    PPG_df["Heure"] = PPG_df["Heure"]/sfreq
# Récupération de l'heure, en seconde de départ
    heureDebut = annot_df.loc[0, "Heure"]
# Modification de la colonne de temps du fichier .EDF
    PPG_df["Heure"] += heureDebut

    return (PPG_df, sfreq)
