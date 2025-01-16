#!/bin/bash

# Déterminer les répertoires par défaut
directories=("resultsStrongScaling/Deisa" "resultsStrongScaling/NoDeisa")

# Déterminer le nombre de nœuds par défaut
node_count="1"

# Vérifier si l'utilisateur a fourni des arguments
if [ $# -eq 1 ]; then
    # Vérifier si l'argument est Deisa ou NoDeisa
    if [ "$1" == "Deisa" ]; then
        directories=("resultsStrongScaling/Deisa")
    elif [ "$1" == "NoDeisa" ]; then
        directories=("resultsStrongScaling/NoDeisa")
    else
        echo "Usage: $0 [Deisa|NoDeisa] [node_count]"
        exit 1
    fi
elif [ $# -eq 2 ]; then
    # Vérifier si l'argument est Deisa ou NoDeisa
    if [ "$1" == "Deisa" ]; then
        directories=("resultsStrongScaling/Deisa")
    elif [ "$1" == "NoDeisa" ]; then
        directories=("resultsStrongScaling/NoDeisa")
    else
        echo "Usage: $0 [Deisa|NoDeisa] [node_count]"
        exit 1
    fi

    # Vérifier si le nombre de nœuds est un entier positif
    if ! [[ "$2" =~ ^[0-9]+$ ]]; then
        echo "Usage: $0 [Deisa|NoDeisa] [node_count]"
        exit 1
    fi

    node_count="$2"
fi

# Parcourir les répertoires par défaut
for directory in "${directories[@]}"; do
    # Obtention du nom du répertoire (Deisa ou NoDeisa)
    dirname=$(basename "$directory")
    
    # Parcourir les sous-répertoires 1, 4 et 16
    for subdir in "$directory"/*; do
        # Obtention du nombre dans le nom du sous-répertoire
        subdirname=$(basename "$subdir")
        num=$(echo "$subdirname" | tr -dc '0-9')
        
        # Vérifier si le nombre de nœuds correspond au nombre spécifié
        if [ "$num" == "$node_count" ]; then
            # Parcourir les fichiers res128.out, res192.out, res256.out et res320.out
            for file in "$subdir"/*; do
                # Obtention de la taille dans le nom du fichier
                filename=$(basename "$file")
                size=$(echo "$filename" | grep -oE '[0-9]+')
                
                # Filtrer les lignes commençant par [RESULT]
                result=$(grep '^\[RESULT\]' "$file")
                
                # Vérifier si le résultat n'est pas vide avant d'afficher
                if [ -n "$result" ]; then
                    echo -e "[$dirname|$num|$size]"
		    echo -e "$result\n"
                fi
            done
        fi
    done
done

