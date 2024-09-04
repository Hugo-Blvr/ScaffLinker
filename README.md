# ScaffLinker v1.0.0

ScaffLinker est un outil de scaffolding qui permet d'assembler des contigs de génomes fragmentés en les organisant en scaffolds, augmentant ainsi la continuité et la qualité des assemblages génomiques. En l'absence de génome de référence ou d'autres données génomiques, ScaffLinker utilise uniquement les assemblages disponibles pour améliorer la structure globale des génomes, essentiel pour les analyses de génomique comparative et de pangénomique.

## Table des matières

- [Dépendances](#dépendances)
- [Installation](#installation)
- [Utilisation](#utilisation)
- [Paramètres](#paramètres)
- [Licence](#licence)


## Dépendances

Pour utiliser ScaffLinker, assurez-vous que les dépendances suivantes sont installées :

- **Système d'exploitation** : Linux, macOS
- **Python 3.x** : ScaffLinker nécessite Python3 pour exécuter le script principal.
- **Pandas** : 
    ```bash
    pip install pandas
    ```
- **Biopython** : 
    ```bash
    pip install biopython
    ```
- **Minimap2** : Doit être installé séparément, téléchargez-le à partir de [GitHub](https://github.com/lh3/minimap2) et suivez les instructions d'installation.


## Installation

Pour installer ScaffLinker, suivez ces étapes :

1. Clonez le dépôt pour obtenir les fichiers nécessaires à l'exécution de ScaffLinker :

    ```bash
    git clone https://github.com/Hugo-Blvr/ScaffLinker.git
    ```

2. Accédez au répertoire du projet :

    ```bash
    cd ScaffLinker
    ```

## Utilisation

### Commande de base

Pour utiliser ScaffLinker, vous avez uniquement besoin d'un dossier contenant des fichiers FASTA assemblés.

Utilisez la commande suivante :

```bash 
./ScaffLinker.sh --dir_assemblies <chemin_vers_votre_répertoire_fasta> [options]
```

Remplacez `<chemin_vers_votre_répertoire_fasta>` par le chemin vers votre répertoire contenant les fichiers FASTA assemblés. Assurez-vous que chaque fichier FASTA a un identifiant unique en préfixe suivi d'un underscore (`_`).

Des répertoires d'assemblages de test sont fournis pour apprendre à utiliser l'outil.

### Exemples d'utilisation test

1. **Utiliser les paramètres par défaut avec le dossier de test fourni** :

   ```bash 
   ./ScaffLinker.sh -d 01_filtered_assemblies
   ```

2. **Spécifier un fichier FASTA unique à échafauder et utiliser plusieurs cœurs pour le processus de mapping** :

   ```bash
   ./ScaffLinker.sh -d 01_filtered_assemblies/ -sf 01_filtered_assemblies/Gd293_filtered50K.fasta -t 5
   ```

3. **Définir un nombre minimum de correspondances et un seuil d'identité de séquence** :

   ```bash
   ./ScaffLinker.sh -d 01_filtered_assemblies --NbMatch 10000 --IdSeq 0.95
   ```

4. **Utiliser des assemblages avec les régions répétées masqués pour le mapping avec plusieurs cœurs** :

   ```bash
   ./ScaffLinker.sh -d 01_filtered_assemblies -dm 01_filtered_assemblies_masked -t 4
   ```

5. **Exécuter ScaffLinker avec un fichier FASTA unique et un nombre de matchs minimum différent** :

   ```bash
   ./ScaffLinker.sh -d 01_filtered_assemblies/ -sf 01_filtered_assemblies/Gd45_filtered50K.fasta -n 8000
   ```


## Paramètres

- `-d`, `--dir_assemblies` : Répertoire contenant les assemblages (obligatoire).
- `-dm`, `--dir_assemblies_masked` : Répertoire contenant les assemblages avec les régions répétées masqués pour le mapping.
- `-n`, `--NbMatch` : Nombre minimum de matchs pour le scaffolding (défaut : 5000).
- `-i`, `--IdSeq` : Seuil d'identité de séquence pour le scaffolding (défaut : 0.90).
- `-o`, `--dir_out` : Répertoire de sortie pour les résultats du scaffolding (défaut : 03_scaffolding).
- `-t`, `--thread` : Nombre de cœurs à utiliser pour le processus de mapping (défaut : 1).
- `-sf`, `--single_fasta` : Spécifie un seul fichier FASTA à échafauder (par défaut : tous les fichiers .fasta).
- `-h`, `--help` : Affiche le message d'aide et quitte.

## Licence
Ce projet est sous licence GNU General Public License v2.0. Consultez le fichier [LICENSE](LICENSE) pour plus de détails.

