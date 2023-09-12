import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from numpy import matlib
import pywt

### Filtrage des artefacts liés au capteur
def filter_rms (PPG_df, seuil, sfreq):
    width = int(sfreq * 0.1) ## Largeur de la fenêtre mobile

    data = (PPG_df.PPG ** 2).rolling(width).sum()
    rms = np.sqrt(data / width)
    masque = np.where(rms < seuil)[0]
    PPG_df_rms = PPG_df.drop(masque)
    return PPG_df_rms

### Detrend - Application d'un filtre passe-haut par transformée en ondelettes
def filtre_PasseHaut_wv (PPG_df):
    wavelet_type='db6' # Type de wavelet
    mode = 'symmetric'

    data = PPG_df.PPG
    DWTcoeffs = pywt.wavedec(data,wavelet_type,mode, level=9, axis=-1)
    DWTcoeffs[0] = np.zeros_like(DWTcoeffs[0])
    filtered_data_haut=pywt.waverec(DWTcoeffs,wavelet_type,mode,axis=-1)
    if len(filtered_data_haut) > len(PPG_df) :
        filtered_data_haut = filtered_data_haut[:-1]
    PPG_df = pd.DataFrame({"Heure" : PPG_df.Heure, "PPG" : filtered_data_haut})

    return PPG_df

### PRE TRAITEMENT DU SIGNAL DE PPG
def PPG_preProcessing (PPG_df, sfreq):
    PPG_df = filter_rms(PPG_df, 5, sfreq)
    PPG_df = filtre_PasseHaut_wv(PPG_df)
    PPG_df.PPG = signal.savgol_filter(PPG_df.PPG, int(sfreq * 0.2), 2) # Fenêtre de 200 ms

     ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(PPG_df)-1, num = len(PPG_df), dtype = int)
    PPG_df = PPG_df.assign(index = index)
    PPG_df.set_index(index, drop = True, inplace = True)
    PPG_df.drop(columns = "index", inplace = True)

    return PPG_df

###

def get_extrema(PPG_df):
    ## Repérage des peaks et des nadirs
    (idx_peaks, properties_peaks) = signal.find_peaks(PPG_df.PPG, height = 0)
    peaks_df = pd.DataFrame({"Heure" : PPG_df.loc[idx_peaks, "Heure"], "PPG" : PPG_df.loc[idx_peaks, "PPG"]})
    PPG_inverted = PPG_df.PPG * (-1)
    (idx_nadir, properties_nadir) = signal.find_peaks(PPG_inverted, height = PPG_inverted.mean())
    nadir_df = pd.DataFrame({"Heure" : PPG_df.loc[idx_nadir, "Heure"], "PPG" : PPG_df.loc[idx_nadir, "PPG"]})
    ## Ajout d'une colonne "peak" / "nadir" sur chaque dataframe
    typeP = np.matlib.repmat(1, len(peaks_df), 1) ## Peaks = + 1
    typeN = np.matlib.repmat(-1, len(nadir_df), 1) ## Nadirs = - 1
    peaks_df = peaks_df.assign(type = typeP)
    nadir_df = nadir_df.assign(type = typeN)
    ## Merge entre peaks_df et nadir_df
    extrema_df = peaks_df.merge(nadir_df, how = "outer", on = ("Heure", "type", "PPG"), sort = True)
    ## Filtrage
        ## Si un peak est suivi d'un peak, on supprime le dernier peak
        ## Si un nadir est suivi d'un nadir, on supprime le dernier nadir
    masque = np.array([])
    for n in range(0,len(extrema_df)-1):
        if (extrema_df.loc[n, "type"]*extrema_df.loc[n+1, "type"]) > 0 :
            masque = np.append(masque, n+1)
    filtered_extrema_df = extrema_df.drop(masque)

    ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(filtered_extrema_df)-1, num = len(filtered_extrema_df), dtype = int)
    filtered_extrema_df = filtered_extrema_df.assign(index = index)
    filtered_extrema_df.set_index(index, drop = True, inplace = True)
    filtered_extrema_df.drop(columns = "index", inplace = True)

    return filtered_extrema_df

### CALCUL PWA
def get_PWA (PPG_df):
    extrema_df = get_extrema(PPG_df)

    PWAheure_np = np.array([])
    PWA_np = np.array([])
    for n in range(0, len(extrema_df)-1):
        if (extrema_df.loc[n, "type"] == 1):
            PWAheure_np = np.append(PWAheure_np, extrema_df.loc[n,"Heure"])
            PWA_np = np.append(PWA_np, extrema_df.loc[n, "PPG"] - extrema_df.loc[n+1, "PPG"])

    PWA_df = pd.DataFrame({'Heure' : PWAheure_np, 'PWA' : PWA_np})
    return PWA_df

###

def filter_PWAnonPhysio (PWA_df):
    ## Distance entre deux cycles cardiaque non physiologique < 0,24 secondes (i.e. HB > 250 / minutes)
    masque = np.array([])
    for n in range(0,len(PWA_df)-1):
        if ((PWA_df.loc[n+1, "Heure"]-PWA_df.loc[n,"Heure"]) < 0.24):
            masque = np.append(masque, n)
    PWA_filtered_df = PWA_df.drop(masque)

    ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(PWA_filtered_df)-1, num = len(PWA_filtered_df), dtype = int)
    PWA_filtered_df = PWA_filtered_df.assign(index = index)
    PWA_filtered_df.set_index(index, drop = True, inplace = True)
    PWA_filtered_df.drop(columns = "index", inplace = True)
    return PWA_df

def filter_PWA_outliers (PWA_df): ## Reprise du code de Monica Betta et al.
    delta_c = abs(np.diff(PWA_df.PWA))
    delta_pp = delta_c[1:]
    delta_c = delta_c[0:-1]
    THR = 3 * np.std(abs(delta_c))
    idx_out = np.where((delta_c > THR)&(delta_pp > THR))[0] + 1
    PWA_df = PWA_df.drop(idx_out)
    return PWA_df

### FILTRAGE DU SIGNAL DE PWA
def filter_PWA (PWA_df):
    PWA_df = filter_PWAnonPhysio (PWA_df)
    PWA_df = filter_PWA_outliers (PWA_df)
    ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(PWA_df)-1, num = len(PWA_df), dtype = int)
    PWA_df = PWA_df.assign(index = index)
    PWA_df.set_index(index, drop = True, inplace = True)
    PWA_df.drop(columns = "index", inplace = True)
    return PWA_df

###
### PRE - TRAITEMENT DU SIGNAL DE PWA
def PWA_Processing (PWA_df):
    PWA_df = pd.DataFrame({"Heure" : PWA_df.Heure, "PWA" : PWA_df.PWA.rolling(int(5)).mean()})
    PWA_df = PWA_df.drop(np.where(np.isnan(PWA_df.PWA))[0])
    PWA_df.PWA = signal.savgol_filter(PWA_df.PWA, 5, 2)
    ## Nettoyage des indices du DataFrame
    index = np.linspace(0, len(PWA_df)-1, num = len(PWA_df), dtype = int)
    PWA_df = PWA_df.assign(index = index)
    PWA_df.set_index(index, drop = True, inplace = True)
    PWA_df.drop(columns = "index", inplace = True)
    return PWA_df

###

def get_PWA_final (PPG_df, sfreq):
    PPG_pp = PPG_preProcessing (PPG_df, sfreq)
    PWA_df = get_PWA (PPG_pp)
    PWA_df = filter_PWA (PWA_df)
    PWA_final_df = PWA_Processing (PWA_df)
    return PWA_final_df

def get_dcc(PWA_df):
    dcc = np.array([])
    for n in range(0, len(PWA_df)-1):
        dcc = np.append(dcc, PWA_df.loc[n+1, "Heure"] - PWA_df.loc[n, "Heure"])

    dcc_df = pd.DataFrame({"Heure": PWA_df.Heure[:-1], "dcc" : dcc})
    return dcc_df
