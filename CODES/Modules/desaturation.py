## Repérage des évènements de désaturation - resaturation à partir du signal de saturation en oxygène dans le sang

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def f_desaturation (sat_df, freq, annot_df, seuil_drop):
    ### PARAMÈTRES -- DÉFINITION D'UN ÉVÈNEMENT DE DÉSATURATION / RESATURATION ###
    seuil_delta_desat = seuil_drop
    seuil_delta_resat = 2/3 # (* delta réel de la désaturation)
    duree_max_desat = 30
    duree_min_event = 10
    duree_max_event = 120

    pas = 1 # Pas de déplacement de la fenêtre temporelle

    desat_df = pd.DataFrame()
    i = 0
    t = 0

    indice_baisse = np.where(np.diff(sat_df.SAT) < 0)[0]

    while ((t < len(indice_baisse))):
        #&(indice_baisse[t] <= len(sat_df)-duree_max_desat*freq)
        n = indice_baisse[t]
        borne_inf = n
        borne_sup = (n + (duree_max_desat*freq))
        if (borne_sup > float(sat_df.index[-1])):
            borne_sup = sat_df.index[-1]
        sat_df_slice = sat_df.SAT[borne_inf : borne_sup]
        ## Calculs du maximum et du minimum sur l'intervalle.
        max_desat = float(sat_df.loc[n, "SAT"])
        heure_max_desat = float(sat_df.loc[n, "Heure"])
        id_max_desat = n
        min_desat = float(sat_df_slice.min(skipna=True))
        id_min = (sat_df_slice.idxmin(skipna=True))
        heure_min = sat_df.loc[id_min, "Heure"]

        if (np.isnan(min_desat) | np.isnan(max_desat) | (heure_max_desat >= heure_min) ):
            t += 1
        else :
            desat = round(max_desat - min_desat) ## arrondi au dessus
            if (desat < seuil_delta_desat):
                t += 1
            else :
                borne_inf = id_min
                borne_sup = (n + (duree_max_event*freq))

                if (borne_sup > float(sat_df.index[-1])):
                    borne_sup = sat_df.index[-1]

                id_slice_min = np.searchsorted(indice_baisse, borne_inf, 'right')
                id_slice_max = np.searchsorted(indice_baisse, borne_sup, 'left')
                indice_baisse_slice = indice_baisse[id_slice_min : id_slice_max]
                k = 0
                flag = 0
                while ((k < len(indice_baisse_slice)) & (flag == 0)) :
                    id_max_resat = indice_baisse_slice[k]
                    #if (isnan(id_max_resat) == 1):
                        #break
                    max_resat = sat_df.loc[id_max_resat, "SAT"]
                    heure_max_resat = sat_df.loc[id_max_resat, "Heure"]
                    resat = max_resat - min_desat
                    duree  = heure_max_resat - heure_max_desat
                    if ((resat >= (seuil_delta_resat)*desat) &(duree_min_event < duree < duree_max_event)):
                        flag = 1
                    k += 1
                if (flag == 0):
                    t += 1
                else :
                    ## Indice pour récupérer le stade de sommeil
                    idx_debut = np.searchsorted(annot_df.Heure, heure_max_desat, "left")
                    if (idx_debut >= len(annot_df)):
                        idx_debut = len(annot_df)-1
                    idx_fin = np.searchsorted(annot_df.Heure, heure_max_resat, "left")
                    if (idx_fin >= len(annot_df)):
                        idx_fin = len(annot_df)-1
                    #if (stade_debut == "Veille")&(stade_fin == "Veille"):
                    #    t += 1
                    #else :
                    ## Mise en place d'un dataframe
                    desat_df.loc[i, "MAX_desat"] = max_desat
                    desat_df.loc[i, "MIN_desat"] = min_desat
                    desat_df.loc[i, "MAX_resat"] = max_resat
                    desat_df.loc[i, "HEURE_debut"] = heure_max_desat # On considère début desat = max desat
                    desat_df.loc[i, "HEURE_min"] = heure_min
                    desat_df.loc[i, "HEURE_fin"] = heure_max_resat
                    desat_df.loc[i, "DUREE_event"] = duree
                    desat_df.loc[i, "DUREE_desat"] = heure_min - heure_max_desat
                    desat_df.loc[i, "DUREE_resat"] = heure_max_resat - heure_min
                    desat_df.loc[i, "DELTA_desat"] = desat
                    desat_df.loc[i, "DELTA_resat"] = resat
                ## Calculs de pente
                    desat_df.loc[i, "PENTE_desat"] = desat /(heure_max_desat - heure_min)
                    desat_df.loc[i, "PENTE_resat"] = resat /(heure_max_resat - heure_min)
                    i += 1
                    t = np.searchsorted(indice_baisse, id_max_resat, 'left')
    return desat_df

def get_trajectoire_desat (desat_df, sat_df):
    i = 0
    m = 0
    heure = {}
    idDESAT = {}
    SAT = {}
    desat_trajectoire_df = pd.DataFrame(columns = ["Heure", "idDESAT", "SAT"])
    for n in range(0, len(desat_df)):
        idxdebut = sat_df[sat_df.Heure == desat_df.loc[n, "HEURE_debut"]].index[0]
        idxfin = (sat_df[sat_df.Heure == desat_df.loc[n, "HEURE_fin"]].index[-1]) + 1
        sat_slice = np.array(sat_df.SAT[idxdebut : idxfin])
        heure_slice = np.array(sat_df.Heure[idxdebut : idxfin])
        idDESAT_slice = np.repeat(i, idxfin - idxdebut)
        new_row = pd.DataFrame({"Heure" : heure_slice, "idDESAT" : idDESAT_slice, "SAT" : sat_slice})
        desat_trajectoire_df = desat_trajectoire_df.append(new_row, ignore_index=True)
        i += 1
    return desat_trajectoire_df
