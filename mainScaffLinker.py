import pandas as pd
import os

pd.set_option('display.max_rows', None)



def Merge_and_filtre(paf_dir, NbMatch, IdSeq):
    """
    Traite les fichiers PAF d'un répertoire donné et filtre les données selon les critères fournis.
    
    Args:
        paf_dir (str): Chemin du répertoire contenant les fichiers PAF.
        NbMatch (int): Seuil minimal pour le nombre de correspondances.
        IdSeq (float): Seuil minimal pour le ratio d'identité de séquence.
    
    Returns:
        tuple: Deux DataFrames pandas, le premier avec les données filtrées selon les critères,
               le second avec les données restantes.
    """
    # Liste des fichiers PAF dans le répertoire
    paf_files = [f"{paf_dir}/{paf}" for paf in os.listdir(paf_dir)]
    # Liste pour stocker les DataFrames individuels
    dataframes = []
    
    # Boucle sur chaque fichier PAF pour le traitement
    for paf_file in paf_files:
        # Chargement du fichier PAF en DataFrame
        df_paf = pd.read_csv(paf_file, sep='\t', header=None).iloc[:, :12]
        # Filtrage sur la qualité de mapping
        df_paf = df_paf[df_paf.iloc[:, 11] > 40].iloc[:, :11]
        # Extraction du nom de l'échantillon à partir du nom de fichier
        paf_name = paf_file.split('/')[-1].split('.')[0].split('_')[-1]
        # Ajout du nom de l'échantillon à la première colonne
        df_paf[0] = f"{paf_name}$" + df_paf[0]
        # Ajout du DataFrame à la liste
        dataframes.append(df_paf)
        
    # Concaténation de tous les DataFrames en un seul
    df = pd.concat(dataframes, ignore_index=True)
    # Calcul de la nouvelle colonne 'IdSeq' (= NbMatch / NbBase)
    df['IdSeq'] = df[9] / df[10]
    # Renommage des colonnes pour plus de clarté
    df.columns = ["Qname", "Qlen", "Qstart", "Qstop", "Strand", "Tname", "Tlen", "Tstart", "Tstop", "NbMatch", "NbBase", "IdSeq"]
    
    # Filtrage des DataFrames selon les critères NbMatch et IdSeq
    df_filtre = df[(df["NbMatch"] >= NbMatch) & (df["IdSeq"] >= IdSeq)]
    df_reste = df[(df["NbMatch"] < NbMatch) | (df["IdSeq"] < IdSeq)]
    
    # Retourner les DataFrames filtrés et non filtrés
    return df_reste.reset_index(drop=True), df_filtre.reset_index(drop=True)



def Ancrage(data):
    """
    Trouve les associations de groupes entre les valeurs de 'Tname' et 'Qname' dans le DataFrame.
    
    Args:
        data (pd.DataFrame): DataFrame contenant les colonnes 'Tname' et 'Qname'.
        
    Returns:
        list: Liste de listes, où chaque sous-liste contient des groupes associés de 'Tname'.
    """

    # Obtenir les valeurs uniques de 'Tname' dans un ensemble pour une recherche plus rapide
    list_tname = set(data['Tname'].unique())
    # Liste pour stocker les groupes associés
    liste_asso = []
    
    # Tant qu'il reste des éléments dans list_tname
    while list_tname:
        # Initialiser le groupe actuel de 'Tname' avec un des 'Tname' restants
        tn = {list_tname.pop()}  # Utilisation de set pour éviter les doublons et recherche rapide
        qn_asos = set(data[data['Tname'].isin(tn)]['Qname'])  # Initialiser les 'Qname' associés
        
        # Boucle pour trouver tous les 'Tname' associés aux 'Qname' actuels
        while True:
            # Trouver tous les 'Tname' associés aux 'Qname' actuels
            tn_asos = set(data[data['Qname'].isin(qn_asos)]['Tname'].unique())
            # Filtrer les 'Tname' déjà traités pour éviter les doublons
            tn_asos -= tn  # Enlever les 'Tname' déjà ajoutés
            
            # Si aucun nouveau 'Tname' n'est trouvé, on arrête la boucle
            if not tn_asos: break

            # Ajouter les nouveaux 'Tname' au groupe actuel
            tn.update(tn_asos)
            # Mettre à jour les 'Qname' associés pour l'itération suivante
            qn_asos = set(data[data['Tname'].isin(tn_asos)]['Qname'])
               
        # Ajouter le groupe trouvé à la liste des associations
        liste_asso.append(list(tn))
        # Mettre à jour list_tname pour retirer les 'Tname' traités
        list_tname -= tn
    
    return liste_asso



def Run(paf_dir, NbMatch, IdSeq):  
    df, df_filtre=Merge_and_filtre(paf_dir, NbMatch, IdSeq)
    associations = Ancrage(df_filtre)
    
    for asso in associations:
        print(asso)


if __name__ == "__main__":
    paf_dir = "02_masked_paf_files/02_paf_files_Gd45"
    NbMatch = 10000
    IdSeq = 0.90
    Run(paf_dir, NbMatch, IdSeq)

    print('Terminée')