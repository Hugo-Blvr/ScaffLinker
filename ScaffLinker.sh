#!/bin/bash

# --help
usage() {
    printf "\nUsage: %s --dir_assemblies <directory> [options]\n\n" "$(basename "$0")"
    printf "Description:\n"
    printf "  ScaffLinker effectue le scaffolding de fichiers FASTA en utilisant un ensemble\n"
    printf "  de paramètres définis par l'utilisateur.\n"
    printf "IMPORTANT : Chaque fichier .fasta doit comporter un identifiant unique en préfixe, suivi d'un underscore ('_').\n\n"
    printf "Options obligatoires:\n"
    printf "  -d, --dir_assemblies          Répertoire contenant les assemblages (obligatoire)\n\n"
    printf "Options:\n"
    printf "  -dm, --dir_assemblies_masked  Répertoire contenant les assemblages avec les régions répétées masqués à utiliser pour le mapping\n"
    printf "  -n,  --NbMatch                Nombre minimum de matchs pour le scaffolding (défaut : 5000)\n"
    printf "  -i,  --IdSeq                  Seuil d'identité de séquence à utiliser pour le scaffolding (défaut : 0.90)\n"
    printf "  -o,  --dir_out                Répertoire de sortie pour les résultats du scaffolding (défaut : 03_scaffolding)\n"
    printf "  -t,  --thread                 Nombre de cœurs à utiliser pour le processus de mapping (défaut : 1)\n"
    printf "  -sf, --single_fasta           Spécifie un seul fichier FASTA à scaffolder (défaut : tous les fichiers .fasta)\n"
    printf "  -h,  --help                   Afficher ce message d'aide et quitter\n"
    printf "\nExemples d'utilisation:\n"
    printf "  %s --dir_assemblies /chemin/vers/assemblages --NbMatch 10000 --IdSeq 0.95\n" "$(basename "$0")"
    printf "  %s -d /chemin/vers/assemblages -dm /chemin/vers/assemblages_masked -t 4\n" "$(basename "$0")"
    printf "\n"
}


# Initialiser les variables avec les valeurs par défaut
initialize_defaults() {
    DIR_ASSEMBLIES=""
    DIR_ASSEMBLIES_MASKED=""
    NBMATCH=5000
    IDSEQ=0.90
    DIR_OUT="03_scaffolding"
    THREAD=1
    SINGLE_FASTA=""
}

# Analyser les arguments en ligne de commande
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -d|--dir_assemblies)
                DIR_ASSEMBLIES="$2"
                shift 2
                ;;
            -dm|--dir_assemblies_masked)
                DIR_ASSEMBLIES_MASKED="$2"
                shift 2
                ;;
            -n|--NbMatch)
                NBMATCH="$2"
                shift 2
                ;;
            -i|--IdSeq)
                IDSEQ="$2"
                shift 2
                ;;
            -o|--dir_out)
                DIR_OUT="$2"
                shift 2
                ;;
            -t|--thread)
                THREAD="$2"
                shift 2
                ;;
            -sf|--single_fasta)
                SINGLE_FASTA="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                printf "Option inconnue: %s\n" "$1" >&2
                usage
                exit 1
                ;;
        esac
    done
}

# Valider les paramètres obligatoires et les options fournies
validate_parameters() {
    if [[ -z "$DIR_ASSEMBLIES" ]]; then
        printf "Erreur: --dir_assemblies est obligatoire.\n" >&2
        usage
        exit 1
    fi

    if [[ "$NBMATCH" -lt 1000 ]]; then
        printf "Erreur: NbMatch doit être au moins 1000.\n" >&2
        exit 1
    fi

    if [[ -n "$SINGLE_FASTA" && ! -f "$SINGLE_FASTA" ]]; then
        printf "Erreur: Le fichier FASTA spécifié n'existe pas: %s\n" "$SINGLE_FASTA" >&2
        exit 1
    fi
}

# Obtenir la liste des fichiers FASTA dans un répertoire donné
retrieve_fasta_files() {
    local dir="$1"
    local fasta_files

    if ! fasta_files=$(find "$dir" -type f -name "*.fasta" 2>/dev/null); then
        printf "Erreur: Aucun fichier FASTA trouvé dans %s.\n" "$dir" >&2
        return 1
    fi

    printf "%s\n" "$fasta_files"
}

# Effectuer le processus de mapping et de scaffolding
perform_mapping_and_scaffolding() {
    local fasta_files
    local target_directory

    # Déterminer le répertoire à utiliser pour le mapping
    if [[ -n "$DIR_ASSEMBLIES_MASKED" && -d "$DIR_ASSEMBLIES_MASKED" ]]; then
        target_directory="$DIR_ASSEMBLIES_MASKED"
    else
        target_directory="$DIR_ASSEMBLIES"
    fi

    # Récupérer les fichiers FASTA à partir du répertoire cible
    if ! fasta_files=$(retrieve_fasta_files "$target_directory"); then
        exit 1
    fi

    # Créer le répertoire de sortie et nettoyer les fichiers existants
    mkdir -p "$DIR_OUT" || { printf "Erreur: Impossible de créer le répertoire de sortie %s.\n" "$DIR_OUT" >&2; exit 1; }
    rm -f "$DIR_OUT"/*

    # Traiter chaque fichier FASTA pour le mapping
    if [[ -n "$SINGLE_FASTA" ]]; then
        local target_file="$SINGLE_FASTA"
        local target_name
        target_name=$(basename "$target_file" | cut -d'_' -f1)
        process_mapping_for_target "$target_file" "$fasta_files" "$target_name"
    else
        for target_file in $fasta_files; do
            local target_name
            target_name=$(basename "$target_file" | cut -d'_' -f1)
            process_mapping_for_target "$target_file" "$fasta_files" "$target_name"
        done
    fi
}

# Traiter le mapping pour chaque fichier FASTA cible
process_mapping_for_target() {
    local target_file="$1"
    local fasta_files="$2"
    local target_name="$3"
    local paf_directory="02_masked_paf_files/02_paf_files_${target_name}"

    # Créer le répertoire pour les fichiers PAF et le nettoyer
    mkdir -p "$paf_directory" || { printf "Erreur: Impossible de créer le répertoire %s.\n" "$paf_directory" >&2; return 1; }
    rm -f "$paf_directory"/*

    printf "Mapping de %s\n" "$target_name"
    for query_file in $fasta_files; do
        if [[ "$query_file" != "$target_file" ]]; then
            local query_name
            query_name=$(basename "$query_file" | cut -d'_' -f1)
            printf "\tMapping de %s et %s ... ...\n" "$target_name" "$query_name"
            if ! minimap2 "$target_file" "$query_file" -o "${paf_directory}/${target_name}_${query_name}.paf" -t "$THREAD" > /dev/null 2> "${paf_directory}/${target_name}_${query_name}.err"; then
                printf "Erreur: Mapping échoué entre %s et %s.\n" "$target_name" "$query_name" >&2
                return 1
            fi
            rm -f "${paf_directory}/${target_name}_${query_name}.err"
        fi
    done

    printf "\tDébut du scaffolding pour %s\n" "$target_name"
    if ! python3 mainScaffLinker.py "$DIR_ASSEMBLIES" "$paf_directory" "$NBMATCH" "$IDSEQ" "$DIR_OUT"; then
        printf "Erreur: Scaffolding échoué pour %s.\n" "$target_name" >&2
        return 1
    fi
}

# Fonction principale
main() {
    initialize_defaults
    parse_arguments "$@"
    validate_parameters
    perform_mapping_and_scaffolding
}

main "$@"
