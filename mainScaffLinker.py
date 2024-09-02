import pandas as pd
import os

#pd.set_option('display.max_rows', None)


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


def Del_repeat(recup, seuil):
    """
    Supprime les lignes du DataFrame où les différences dans les coordonnées de début et de fin
    pour Tname ou Qname sont inférieures à un seuil prédéfini, indiquant un 'repeat'.

    Args:
        recup (pd.DataFrame): DataFrame contenant les informations de correspondance.
        seuil (int): Le seuil pour la différence maximale pour déterminer un 'repeat'.

    Returns:
        pd.DataFrame: DataFrame mis à jour sans les 'repeats'.
    """
    df = recup.copy()

    # Grouper par 'Tname' et 'Qname' pour l'analyse
    grouped = df.groupby(['Tname', 'Qname'])

    # Initialiser les listes pour les indices et les contigs à supprimer
    false_contig = set()
    indexes_to_drop = []

    # Itération sur les groupes pour identifier les 'repeats'
    for (tname, qname), group in grouped:
        if len(group) > 1:
            # Calcul des différences pour déterminer les 'repeats'
            tstart_diff = group['Tstart'].max() - group['Tstart'].min()
            tstop_diff = group['Tstop'].max() - group['Tstop'].min()
            qstart_diff = group['Qstart'].max() - group['Qstart'].min()
            qstop_diff = group['Qstop'].max() - group['Qstop'].min()

            target_repeat = (tstart_diff < seuil and tstop_diff < seuil)
            query_repeat = (qstart_diff < seuil and qstop_diff < seuil)

            # Marquer les contigs ou queries comme 'false' si répétition
            if target_repeat or query_repeat:
                if target_repeat: false_contig.add(tname)
                if query_repeat: false_contig.add(qname)
                indexes_to_drop.extend(group.index)

    # Vérifier les groupes uniques et supprimer si marqués comme 'false'
    for (tname, qname), group in grouped:
        if len(group) == 1 and (tname in false_contig or qname in false_contig): 
            indexes_to_drop.extend(group.index)

    # Supprimer les indices identifiés du DataFrame
    if indexes_to_drop: df.drop(indexes_to_drop, inplace=True)

    return df.reset_index(drop=True)


def Reverse(recup, direction):
    """
    Modifie les coordonnées de Tstart, Tstop, Qstart, Qstop, et la direction de 'Strand'
    en fonction de la direction de référence et des contigs inversés.

    Args:
        recup (pd.DataFrame): DataFrame contenant les informations de correspondance.
        direction (list): Liste contenant quatre listes : Tsens_Tref, Tinv_Tref, Qsens_Tref, Qinv_Tref.

    Returns:
        pd.DataFrame: DataFrame mis à jour avec les coordonnées et directions modifiées.
    """
    # Créer une copie du DataFrame pour ne pas modifier l'original
    data = recup.copy()

    # Ajouter les colonnes 'ReverseT' et 'ReverseQ' initialisées à False
    data["ReverseT"] = False
    data["ReverseQ"] = False

    # Décompacter les listes de direction
    Tsens_Tref, Tinv_Tref, Qsens_Tref, Qinv_Tref = direction

    # Masques pour chaque condition
    mask_Tsens_Qinv = data["Tname"].isin(Tsens_Tref) & data["Qname"].isin(Qinv_Tref)
    mask_Tinv_Qsens = data["Tname"].isin(Tinv_Tref) & data["Qname"].isin(Qsens_Tref)
    mask_Tinv_Qinv = data["Tname"].isin(Tinv_Tref) & data["Qname"].isin(Qinv_Tref)

    # Calculer les nouvelles valeurs pour le masque Tsens_Tref et Qinv_Tref
    new_Qstart_Tsens_Qinv = data.loc[mask_Tsens_Qinv, "Qlen"] - data.loc[mask_Tsens_Qinv, "Qstop"]
    new_Qstop_Tsens_Qinv = data.loc[mask_Tsens_Qinv, "Qlen"] - data.loc[mask_Tsens_Qinv, "Qstart"]
    # Appliquer les modifications pour le masque Tsens_Tref et Qinv_Tref
    data.loc[mask_Tsens_Qinv, "Qstart"] = new_Qstart_Tsens_Qinv
    data.loc[mask_Tsens_Qinv, "Qstop"] = new_Qstop_Tsens_Qinv
    data.loc[mask_Tsens_Qinv, "ReverseQ"] = True
    data.loc[mask_Tsens_Qinv, "Strand"] = data.loc[mask_Tsens_Qinv, "Strand"].apply(lambda x: '-' if x == '+' else '+')

    # Calculer les nouvelles valeurs pour le masque Tinv_Tref et Qsens_Tref
    new_Tstart_Tinv_Qsens = data.loc[mask_Tinv_Qsens, "Tlen"] - data.loc[mask_Tinv_Qsens, "Tstop"]
    new_Tstop_Tinv_Qsens = data.loc[mask_Tinv_Qsens, "Tlen"] - data.loc[mask_Tinv_Qsens, "Tstart"]
    # Appliquer les modifications pour le masque Tinv_Tref et Qsens_Tref
    data.loc[mask_Tinv_Qsens, "Tstart"] = new_Tstart_Tinv_Qsens
    data.loc[mask_Tinv_Qsens, "Tstop"] = new_Tstop_Tinv_Qsens
    data.loc[mask_Tinv_Qsens, "ReverseT"] = True
    data.loc[mask_Tinv_Qsens, "Strand"] = data.loc[mask_Tinv_Qsens, "Strand"].apply(lambda x: '-' if x == '+' else '+')

    # Calculer les nouvelles valeurs pour le masque Tinv_Tref et Qinv_Tref
    new_Tstart_Tinv_Qinv = data.loc[mask_Tinv_Qinv, "Tlen"] - data.loc[mask_Tinv_Qinv, "Tstop"]
    new_Tstop_Tinv_Qinv = data.loc[mask_Tinv_Qinv, "Tlen"] - data.loc[mask_Tinv_Qinv, "Tstart"]
    new_Qstart_Tinv_Qinv = data.loc[mask_Tinv_Qinv, "Qlen"] - data.loc[mask_Tinv_Qinv, "Qstop"]
    new_Qstop_Tinv_Qinv = data.loc[mask_Tinv_Qinv, "Qlen"] - data.loc[mask_Tinv_Qinv, "Qstart"]
    # Appliquer les modifications pour le masque Tinv_Tref et Qinv_Tref
    data.loc[mask_Tinv_Qinv, "Tstart"] = new_Tstart_Tinv_Qinv
    data.loc[mask_Tinv_Qinv, "Tstop"] = new_Tstop_Tinv_Qinv
    data.loc[mask_Tinv_Qinv, "Qstart"] = new_Qstart_Tinv_Qinv
    data.loc[mask_Tinv_Qinv, "Qstop"] = new_Qstop_Tinv_Qinv
    data.loc[mask_Tinv_Qinv, "ReverseT"] = True
    data.loc[mask_Tinv_Qinv, "ReverseQ"] = True

    return data


def Verification(reverse):
    """
    Agrège les informations par 'Tname' et 'Qname', calcule les couvertures et filtre les résultats
    selon des critères de couverture et de positions de début et de fin.

    Args:
        reverse (pd.DataFrame): DataFrame contenant les correspondances inversées.

    Returns:
        pd.DataFrame: DataFrame agrégé et filtré selon les critères spécifiés.
    """
    # Créer une copie du DataFrame
    rv = reverse.copy()
    
    # Agréger les données par 'Tname' et 'Qname'
    grouped = rv.groupby(['Tname', 'Qname']).agg({
        'Qlen': 'mean',
        'Qstart': 'min',
        'Qstop': 'max',
        'Tlen': 'mean',
        'Tstart': 'min',
        'Tstop': 'max',
        'NbMatch': 'sum',
        'IdSeq': 'mean',
        'ReverseT': 'first',
        'ReverseQ': 'first'
    }).reset_index()

    # Calculer 'Qcover' et 'Tcover'
    grouped['Qcover'] = grouped['NbMatch'] / (grouped['Qstop'] - grouped['Qstart'])
    grouped['Tcover'] = grouped['NbMatch'] / (grouped['Tstop'] - grouped['Tstart'])

    # Filtrer selon 'Qcover' et 'Tcover' > 0.3
    filtered = grouped[(grouped['Qcover'] > 0.3) & (grouped['Tcover'] > 0.3)]

    # Filtrer selon les critères de positions
    filtered = filtered[((filtered['Qstart'] < 100000) | (filtered['Qlen'] - filtered['Qstop'] < 100000)) &
                        ((filtered['Tstart'] < 100000) | (filtered['Tlen'] - filtered['Tstop'] < 100000))]

    return filtered


def Run(paf_dir, NbMatch, IdSeq, display = True):  
    df, df_filtre=Merge_and_filtre(paf_dir, NbMatch, IdSeq)
    associations = Ancrage(df_filtre)
    
    for asso in associations:
            
        data_asso = df_filtre[df_filtre['Tname'].isin(asso)].reset_index(drop=True)
        recup1 = Recup_match(df, data_asso,IdSeq-0.15) 
        direction = Direction_assignment(recup1)
        recup2 = Recup_match(df, data_asso,0.6)
        
        recup_less_del = Del_repeat(recup2, 3000)
        new_ancrage = Ancrage(recup_less_del)

        if not new_ancrage: continue
        if len(new_ancrage) > 1 :
            recup_less_del = recup_less_del[recup_less_del['Tname'].isin(new_ancrage[0])]
            for i in new_ancrage[1:]: 
                associations.append(i)
        
        reverse = Reverse(recup_less_del, direction).sort_values(by = ['Tname','Qname','Tstart']).reset_index(drop = True)
        verif = Verification(reverse)
        
        if verif.empty: continue 
        new_ancrage = Ancrage(verif)
        if len(verif) > 1 :
            verif = verif[verif['Tname'].isin(new_ancrage[0])].reset_index(drop = True)
            for i in new_ancrage[1:]: associations.append(i)   
        
        if display : 
            print('\n', '*'*50,new_ancrage,'\n')
            print('\n',reverse,'\n',verif)



if __name__ == "__main__":
    paf_dir = "02_masked_paf_files/02_paf_files_Gd45"
    NbMatch = 10000
    IdSeq = 0.90
    Run(paf_dir, NbMatch, IdSeq)

    print('Terminée')