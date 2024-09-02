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


def Recup_match(df, df_asso, id_seq):
    """
    Récupère les lignes correspondant aux associations spécifiques entre 'Qname' et 'Tname'
    à partir de 'df' et 'df_asso' et applique un filtre basé sur 'NbMatch' et 'IdSeq'.
    
    Args:
        _df (pd.DataFrame): DataFrame contenant les données générales de correspondance.
        df_asso (pd.DataFrame): DataFrame contenant les associations spécifiques de 'Qname' et 'Tname'.
        id_seq (float): Seuil minimal pour le ratio d'identité de séquence 'IdSeq'.
    
    Returns:
        pd.DataFrame: DataFrame filtré contenant les lignes correspondant aux associations spécifiques.
    """

    # Filtrer 'df' pour ne conserver que les lignes avec les mêmes 'Qname' et 'Tname' que dans 'df_asso'
    filtered_df = df[df[['Qname', 'Tname']].apply(tuple, axis=1).isin(df_asso[['Qname', 'Tname']].apply(tuple, axis=1))]
    # Concaténer les lignes filtrées de 'df' à 'df_asso'
    recup_match = pd.concat([df_asso, filtered_df], ignore_index=True)
    # Appliquer les filtres sur 'NbMatch' et 'IdSeq'
    recup_match = recup_match[(recup_match["NbMatch"] > 1000) & (recup_match["IdSeq"] > id_seq)]
    
    # Réinitialiser l'index avant de retourner le DataFrame
    return recup_match.reset_index(drop=True)


def Direction_assignment(data_asso):
    """
    Assigne les directions de correspondance entre contigs 'Tname' et 'Qname' basées sur le nombre de correspondances
    et leur orientation. Le contig avec le plus grand nombre de correspondances est utilisé comme référence.
    
    Args:
        data_asso (pd.DataFrame): DataFrame contenant les colonnes 'Tname', 'Qname', 'NbMatch', et 'Strand'.
    
    Returns:
        list: Listes des contigs avec leur orientation relative au contig de référence.
    """
    
    # Calculer le nombre total de correspondances par 'Tname', 'Qname' et 'Strand'
    coverage_data = (
        data_asso
        .groupby(['Tname', 'Qname', 'Strand'], as_index=False)
        .agg({'NbMatch': 'sum'})
    )
    # Trouver le strand avec le plus de correspondances pour chaque paire ('Tname', 'Qname')
    strand_max = (
        coverage_data
        .loc[coverage_data.groupby(['Tname', 'Qname'])['NbMatch'].idxmax()]
        .reset_index(drop=True)
    )
    
    # Identifier le contig avec le plus grand nombre de correspondances
    T_max_Nbmatch = strand_max.groupby('Tname')['NbMatch'].sum().idxmax()
    
    # Initialiser les listes de contigs et leurs orientations
    Tsens_Tref = [T_max_Nbmatch]
    Tinv_Tref, Qsens_Tref, Qinv_Tref = [], [], []
    
    # Créer un ensemble unique de tous les contigs 'Tname' et 'Qname' présents dans 'strand_max'
    unique_contig = set(strand_max['Tname']).union(set(strand_max['Qname']))
    # Retirer 'T_max_Nbmatch' de l'ensemble 'unique_contig' pour éviter de le traiter dans les boucles d'exploration suivantes
    unique_contig.discard(T_max_Nbmatch)
    # Liste pour stocker les contigs à explorer
    listeT_explore = [T_max_Nbmatch]
    
    # Boucle principale pour explorer les contigs et assigner les directions
    while unique_contig:   
        # Liste pour stocker les 'Qname' qui seront explorés dans cette itération
        listeQ_explore = []
        # Filtrer le DataFrame 'strand_max' pour obtenir les lignes où 'Tname' est dans 'listeT_explore'
        Trow = strand_max[strand_max['Tname'].isin(listeT_explore)]
        # Trouver les lignes avec le nombre maximal de correspondances ('NbMatch') pour chaque 'Qname'
        Tmatch_max = Trow.loc[Trow.groupby('Qname')['NbMatch'].idxmax()]
        
        # Parcourir chaque ligne du DataFrame 'Tmatch_max' pour déterminer les orientations
        for _, row in Tmatch_max.iterrows():
            target, query, strand = row['Tname'], row['Qname'], row['Strand']
            
            # Vérifier si le 'Qname' n'a pas encore été exploré
            if query in unique_contig:
                # Ajouter 'Qname' à la liste des éléments à explorer dans la prochaine itération
                listeQ_explore.append(query)
                
                # Assigner l'orientation de 'Qname' par rapport au contig de référence
                if strand == '+':
                    if target in Tsens_Tref: Qsens_Tref.append(query)
                    else: Qinv_Tref.append(query)
                else:
                    if target in Tsens_Tref: Qinv_Tref.append(query)
                    else: Qsens_Tref.append(query)
        
        # Retirer les 'Tname' explorés de l'ensemble des contigs uniques pour éviter les doublons
        unique_contig.difference_update(listeT_explore)
        # Réinitialiser la liste des 'Tname' à explorer pour la prochaine itération
        listeT_explore = []
        # Filtrer le DataFrame 'strand_max' pour obtenir les lignes où 'Qname' est dans 'listeQ_explore'
        Qrow = strand_max[strand_max['Qname'].isin(listeQ_explore)]
        # Trouver les lignes avec le nombre maximal de correspondances ('NbMatch') pour chaque 'Tname'
        Qmatch_max = Qrow.loc[Qrow.groupby('Tname')['NbMatch'].idxmax()]
        
        # Parcourir chaque ligne du DataFrame 'Qmatch_max' pour déterminer les orientations
        for _, row in Qmatch_max.iterrows():
            target, query, strand = row['Tname'], row['Qname'], row['Strand']
            
            # Vérifier si le 'Tname' n'a pas encore été exploré
            if target in unique_contig:
                # Ajouter 'Tname' à la liste des éléments à explorer dans la prochaine itération
                listeT_explore.append(target)
                
                # Assigner l'orientation de 'Tname' par rapport au contig de référence
                if strand == '+':
                    if query in Qsens_Tref: Tsens_Tref.append(target)
                    else: Tinv_Tref.append(target)
                else:
                    if query in Qsens_Tref: Tinv_Tref.append(target)
                    else: Tsens_Tref.append(target)
        
        # Retirer les 'Qname' explorés de l'ensemble des contigs uniques pour éviter les doublons
        unique_contig.difference_update(listeQ_explore)
    
    # Retourner les listes de contigs classifiés par leurs orientations relatives
    return [Tsens_Tref, Tinv_Tref, Qsens_Tref, Qinv_Tref]


def Run(paf_dir, NbMatch, IdSeq, display = True):  
    df, df_filtre=Merge_and_filtre(paf_dir, NbMatch, IdSeq)
    associations = Ancrage(df_filtre)
    
    for asso in associations:
            
        data_asso = df_filtre[df_filtre['Tname'].isin(asso)].reset_index(drop=True)
        recup1 = Recup_match(df, data_asso,IdSeq-0.15) 
        direction = Direction_assignment(recup1)
            
        if display : 
            print('\n', '*'*50,'\n',asso)
            print('\n',direction)



if __name__ == "__main__":
    paf_dir = "02_masked_paf_files/02_paf_files_Gd45"
    NbMatch = 10000
    IdSeq = 0.90
    Run(paf_dir, NbMatch, IdSeq)

    print('Terminée')