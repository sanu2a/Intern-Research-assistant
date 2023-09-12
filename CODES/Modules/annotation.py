## Importation des données d'annotation relatives aux signaux de PSG

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

# Le fichier d'annotation importé contient les colonnes suivantes :
    ## Heure
    ## Stade
    ## Position
    ## Evènement
    ## Durée
    ## Temps écoulé (pour les dispositifs DELTAMED)
# Les stades, les positions et les durées doivent être inscrites avec le même vocabulaire que celui utilisé dans le scoring DELTAMED (Cf. Importation des fichiers d'annotation Remlogic)

## Fonctions utiles
def heureToSecondRTF(rtf_df): ## Convertir les heures en secondes dans le fichier rtf
    for n in range (0,len(rtf_df)):
        time = dt.datetime.strptime(rtf_df.loc[n, "Heure"], "%Hh%Mm%Ss")
        if time.hour > 12: # Prise en compte du changement de jour
            rtf_df.loc[n,"Heure"] = time.second + time.minute*60 + time.hour*3600
        else:
            rtf_df.loc[n, "Heure"] = time.second + time.minute*60 + (time.hour+24)*3600
    return rtf_df

def heureToSecondTXT(txt_df): ## Convertir les heures en secondes dans le fichier txt
    for n in range (0,len(txt_df)):
        time = dt.datetime.strptime(txt_df.loc[n, "Heure"], "%H:%M:%S")
        if time.hour > 12: # Prise en compte du changement de jour
            txt_df.loc[n, "Heure"] = time.second + time.minute*60 + time.hour*3600
        else:
            txt_df.loc[n, "Heure"] = time.second + time.minute*60 + (time.hour+24)*3600
    return txt_df

def dureeToSecond(rtf_df): ## Convertir les durées en secondes dans le fichier rtf
    for n in range (0,len(rtf_df)):
        if ("00:" in rtf_df.loc[n, "Durée"]): ## TROUVER UNE MÉTHODE PLUS ROBUSTE !
            time = dt.datetime.strptime(rtf_df.loc[n, "Durée"], "%H:%M:%S")
            rtf_df.loc[n, "Durée"] = time.second + time.minute*60 + time.hour*3600
        else :
            rtf_df.loc[n, "Durée"] = np.nan
    return rtf_df

## Importation des annotations DELTAMED

def annotationDeltamed(path, patient):
     # Importation du fichier .TXT
    dataTXT_df = pd.read_table(glob.glob(path + patient + "*.TXT")[0],
                          skiprows=5, sep='\t', encoding='latin',
                          index_col=False, names = ["Heure", "Stade"])
    # Importation du fichier .RTF
    with open(glob.glob(path + patient + "*.rtf")[0], 'r') as file:
        sample_text = file.read()
        text = rtf_to_text(sample_text)
        x = re.sub(r"{\*?\\.+(;})|\s?\\[A-Za-z0-9]+|\s?{\s?\\[A-Za-z0-9]+\s?|\s?}\s?", ";", text)
        lines = [line.split(' ') for line in x.split('\n')[15:-3]]
        res = [[el for el in sub if el != ''] for sub in lines]
        #for row in res:
        #    [row.pop() for i in range(7, len(row))]
          
        dataRTF_df = pd.DataFrame(res, columns=['index', 'Temps_écoulé', 'Heure', 'Durée', 'Evènement', 'compl1', 'compl2'])
        dataRTF_df = dataRTF_df.drop(['index'], axis = 1)

    for n in range(0, len(dataRTF_df)):
        if dataRTF_df.loc[n, "compl2"] is not None :
            dataRTF_df.loc[n, "compl1"] = dataRTF_df.loc[n, "compl1"] + " " + dataRTF_df.loc[n, "compl2"]
        if dataRTF_df.loc[n, "compl1"] is not None :
            dataRTF_df.loc[n, "Evènement"] = dataRTF_df.loc[n, "Evènement"] + " " + dataRTF_df.loc[n, "compl1"]
    dataRTF_df.drop(["compl1", "compl2"], 1, inplace = True)

    # Mise en place d'une colonne de position distinctes de la colonne évènement
    positions_df = pd.DataFrame.copy(dataRTF_df["Evènement"])
    positions_df.columns = ["Position"]
    dataRTF_df.insert(0, "Position", positions_df)

    ## Passage des heures de secondes.
    heureToSecondTXT (dataTXT_df)
    heureToSecondRTF (dataRTF_df)
    dureeToSecond (dataRTF_df)

    ## Jointure des fichiers txt et rtf par la colonne "Heure"
    fullData_df = dataTXT_df.merge(dataRTF_df, on = 'Heure', how = 'outer')

    ## Nettoyage des annotations
    positions = ["Dorsal","Latéral", "Latéral Gauche","Latéral Droit","Assis" ,"Ventral"] ## Annotation de positions
    stades = ["Veille", "Stade 1","Stade 2","Stade 3","S. Paradoxal"] ## Annotation de stades
    events = ["Hypopnée", "hypopnée Centrale","Apnée", "Apnée Centrale", ## Annotation des évènements
    "Apnée Obstructive", "Apnée Mixte","Désaturation",
    "hypoventilation alvéolaire","Respiration périodique","Limitation de débit", "LUMIERE ETEINTE",
    "Arousal d'origine respiratoire", "Mouvement + arousal" ,"PLM droit", "PLM Gauche",
    "Arousal non spécifique", "Respiration", "Ronflements simples"]

    memostade = np.nan
    memopos = np.nan
    for n in range(0, len(fullData_df)):
        if (fullData_df.loc[n, "Stade"] not in stades): ## Completion de la colonne de stades de sommeil
            fullData_df.loc[n, "Stade"] = np.nan
        if (type(fullData_df.loc[n, "Stade"]) == str):
            memostade = fullData_df.loc[n, "Stade"]
        else :
            fullData_df.loc[n, "Stade"] = memostade
        if (fullData_df.loc[n, "Position"] not in positions): ## Completion de la colonne des positions
            fullData_df.loc[n, "Position"] = np.nan
        if (type(fullData_df.loc[n, "Position"]) == str):
            memopos = fullData_df.loc[n, "Position"]
        else :
            fullData_df.loc[n, "Position"] = memopos
        if (fullData_df.loc[n, "Evènement"] not in events):
            fullData_df.loc[n, "Evènement"] = np.nan

    ## Nettoyage du dataframe
    fullData_df.sort_values("Heure", inplace = True)
    fullData_df.drop_duplicates(inplace = True)
    ## Nettoyage des indices du DataFrame
    fullData_df.index = pd.RangeIndex(0,len(fullData_df),1, dtype=np.int32)

    ## Organisation des colonnes
    fullData_df = fullData_df.reindex(columns=['Temps_écoulé','Heure','Stade','Position','Evènement','Durée'])

    ## Optimisation de la mémoire
    fullData_df.Heure = np.float32(fullData_df.Heure)

    return fullData_df

## Importation des annotations REMLOGIC

## Modification des appelations pour utiliser le vocabulaire de DELTAMED
def stade(st):
    switcher = {
    "SLEEP-S0": "Veille",
    "SLEEP-S1": "Stade 1",
    "SLEEP-S2": "Stade 2",
    "SLEEP-S3": "Stade 3",
    "SLEEP-REM": "S. Paradoxal"
    }
    return switcher.get(st)

def position(pos):
    switcher = {
    "POSITION-SUPINE":  "Dorsal",
    "POSITION-LEFT": "Latéral Gauche",
    "POSITION-RIGHT" :  "Latéral Droit",
    "POSITION-UPRIGHT": "Assis" ,
    "POSITION-PRONE": "Ventral",
    }
    return switcher.get(pos)

def event(ev):
    switcher = {
    "HYPOPNEA" : "Hypopnée",
    "HYPOPNEA-CENTRAL" : "hypopnée Centrale",
    "APNEA" : "Apnée",
    "APNEA-CENTRAL" : "Apnée Centrale",
    "APNEA-OBSTRUCTIVE" : "Apnée Obstructive",
    "APNEA-MIXED" : "Apnée Mixte",
    "DESAT" : "Désaturation",
    "HVA" : "hypoventilation alvéolaire",
    "respiration pÈriodique" : "Respiration périodique",
    "LDI": "Limitation de débit",
    "LIGHTS-OFF" : "LUMIERE ETEINTE",
    "AROUSAL-RESP" : "Arousal d'origine respiratoire",
    "AROUSAL-LM" :  "Mouvement + arousal" ,
    "PLM-RRLM" :  "PLM droit",
    "PLM-LM":  "PLM Gauche",
    "AROUSAL-SPONT" :  "Arousal non spécifique",
    "RESP-BREATH-I" :    "Respiration",
    "SNORE-SINGLE": "Ronflements simples"
    }
    return switcher.get(ev)

def annotationRemlogic(path, patient):
    ## Chargement des données
    dataTXT_df = pd.read_table(glob.glob(path + patient + "*.txt")[0],
                        sep='\t', encoding='latin', skiprows = 207,
                        names=['Stade', 'Position', 'Heure', 'Evènement', 'Durée'], header = 0)
    ## Standardisation des appelations
    dataTXT_df.Stade = dataTXT_df.Stade.apply(lambda x : stade(x))
    dataTXT_df.Position = dataTXT_df.Position.apply(lambda x : position(x))
    dataTXT_df.Evènement = dataTXT_df.Evènement.apply(lambda x : event(x))

    ## Completion des stades de sommeil et des positions
    memostade = np.nan
    memopos = np.nan
    for n in range(0, len(dataTXT_df)):
        if (type(dataTXT_df.loc[n, "Stade"]) == str):
            memostade = dataTXT_df.loc[n, "Stade"]
        else :
            dataTXT_df.loc[n, "Stade"] = memostade

        if (type(dataTXT_df.loc[n, "Position"]) == str):
            memopos = dataTXT_df.loc[n, "Position"]
        else :
            dataTXT_df.loc[n, "Position"] = memopos
    dataTXT_df = heureToSecondTXT(dataTXT_df)

    ## Nettoyage du DataFrame
    dataTXT_df.sort_values("Heure", inplace = True)
    dataTXT_df.drop_duplicates(inplace = True)

    ## Nettoyage des indices du DataFrame
    dataTXT_df.index = pd.RangeIndex(0,len(dataTXT_df),1, dtype=np.int32)

    dataTXT_df = dataTXT_df.reindex(columns=['Heure','Stade','Position','Evènement','Durée'])

    ## Optimisation de la mémoire
    dataTXT_df.Heure = np.float32(dataTXT_df.Heure)

    return dataTXT_df

## IMPORTATION DES DONNÉES

## Remarque : les annotations DELTAMED sont composées d'un fichier .TXT et d'un fichier .rtf_df
## les annotations Remlogic sont composées d'un fichier .txt uniquement
def loadAnnotation(path, patient):
    #if ((glob.glob(path + patient + "*.txt") != [])) : ## Si il n'y a pas de fichier .txt
    if any([file.endswith('.txt') for file in os.listdir(path + patient )]):
        file_list = os.listdir(path + patient )
        filename = file_list[[file.endswith('.txt') for file in os.listdir(path + patient )].index(True)]
        #with open(glob.glob(path + patient + "*.txt")[0], 'rb') as file : ## Ouverture fichier .txt
        with open(path+patient+filename, 'rb') as file : ## Ouverture fichier .txt
            sample_text = file.read()
        if (str(sample_text).find("RemLogic") != -1) : # Présence du mot "RemLogic"
            data = annotationRemlogic(path, patient) # Chargement data avec la fonction adaptée au logiciel Remlogic
    else : ## Si il y a un fichier .txt
        file_list = os.listdir(path + patient )
        filename = file_list[[file.endswith('.rtf') for file in os.listdir(path + patient )].index(True)]
        #with open(glob.glob(path + patient + "*.txt")[0], 'rb') as file : ## Ouverture fichier .txt
        with open(path+patient+filename, 'rb') as file : ## Ouverture fichier .txt
        #with open(glob.glob(path + patient + "*.rtf")[0], 'r') as file : # Ouverture fichier .rtf
            sample_text = file.read()
        if (sample_text.find(b"DELTAMED")!= -1) : # Présence du mot "DELTAMED"
            data = annotationDeltamed(path, patient) # Chargement data avec la fonction adaptée au logiciel Remlogic
    return data

## Lors de l'ajout d'un nouveau logiciel il faut :
    ## Créer une fonction permettant l'importation des données :
        ## contenant au moins les 5 colonnes suivantes : ['Heure','Stade','Position','Evènement','Durée']
        ## avec les noms des stades, évènements et positions conforment à l'orthographe utilisé par le logiciel DELTAMED
        ## avec les heures et la durée des évènements en secondes
    ## Appeler cette fonction dans la fonction loadAnnotation avec un moyen robuste pour identifier le logiciel
    ## Effectuer les modifications nécessaires dans le module edf.py

## Fonctions complémentaires
def get_sleepstage (annot_df, wake, rem) :
    Wake = ["Veille"]
    REM = ["S. Paradoxal"]
    NREM = ["Stade 1", "Stade 2", "Stade 3"]
    sleep_np = np.array([])
    for n in annot_df.Stade:
        if n in Wake :
            sleep_np = np.append(sleep_np, wake)
        elif n in NREM :
            sleep_np = np.append(sleep_np, (wake + rem)/2)
        else :
            sleep_np = np.append(sleep_np, rem)
    sleep_df = pd.DataFrame({"Heure" : annot_df.Heure, "Sleep" : sleep_np})
    return sleep_df
