## apnee-git

### Répertoire "biblio"

Le répertoire "biblio" contient les différents articles de recherche utilisés dans l'ensemble du stage ainsi que d'autres documents utiles.

### Répertoire "Codes"

Le répertoire "code" contient le code source Guillemette utilisé.

### Répertoire "code"

Le répertoire "code" contient le code source principal effectué tout au long du stage.

#### Documentation

- [code/analyse.R](code/analyse.R) : Un brouillon utilisé lors des premières analyses faites sur les données au début du stage.

- [code/pretraitement.py](code/pretraitement.py) : Le code utilisé pour extraire les indicateurs, détecter les hyperventilations et reannoter les différents événements sur l'ensemble de la base de données Mars.

- [code/pobm_use.py](code/pobm_use.py) : Le code utilisé pour extraire les désaturations à l'aide de la bibliothèque pobm mentionnée dans l'article de Levy et al (voir répertoire "biblio").

- [code/Multinomreg_hddac.Rmd](code/Multinomreg_hddac.Rmd) : Contient l'ensemble de l'analyse effectuée sur l'ensemble de la base de données (ACP, boxplots, sélection des variables, etc.), ainsi que la modélisation et les modèles de régression multinomiale et les algorithmes HDDA et HDDC, ainsi que les différents résultats des modèles construits.

- [code/spacem3_files.py](code/spacem3_files.py) : Le code utilisé pour créer l'ensemble des fichiers nécessaires pour faire fonctionner le logiciel spacem3.


### Répertoire "Data"

Le répertoire "Data" contient les données utilisées et les données générées durant le stage.

#### Documentation

- [Data/mars](Data/mars) : Ce répertoire contient l'ensemble des données pour chaque patient extraites lors de la polysomnographie (Ces derniers sont supprimées du repertoire car il sont tres volumineux). Un ensemble de fichiers a été ajouté par la suite dans chaque dossier (les données relatives à chaque patient indépendamment, y compris les désaturations, les hyperventilations, le signal de saturation, les annotations, etc.), ainsi que les données relatives à l'ensemble des patients, comme les niveaux de base. Un fichier CSV contenant la totalité des événements et leurs types, ainsi que leurs indicateurs, est également inclus, tout comme les fichiers utilisés par le logiciel spacem3.
- Voici un lien dropbox vers les fihciers necessaires, il faut l'extraire dans Data/mars : https://www.dropbox.com/scl/fi/5yzzzuuovw8x2fwyedskq/mars.zip?rlkey=m3y5b7e3n4vt8214k04naqi6z&dl=0.


### Rapport

Le lien vers le rapport de stage est : https://www.overleaf.com/read/ctykdqjvybrt
