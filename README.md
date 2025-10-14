# TP Ondes – Résolution numérique de l’équation des ondes

Auteure : Maria Saad
ENSTA Paris – Novembre 2024
Langage : MATLAB

## Objectif

Ce projet met en œuvre la méthode des éléments finis (MEF) pour résoudre numériquement l’équation des ondes en régime temporel.
Le travail couvre un problème stationnaire et un schéma temporel explicite (saute-mouton). Un schéma implicite est également fourni.

## Contenu du dépôt

TP_Ondes_MariaSaad/
├─ src/                # scripts MATLAB (.m)
├─ data/               # maillages Gmsh et données d’entrée
├─ docs/               # rapport et énoncé du TP
├─ README.md
├─ .gitignore

## Méthodologie (résumé)

1. Problème stationnaire : éléments finis P1, assemblage M et K, Dirichlet homogène, validation par solution analytique.
2. Schéma explicite : condition CFL, condensation de masse, suivi des énergies.
3. Schéma implicite (optionnel) : résolution via factorisation de Cholesky.

## Utilisation

1. Ouvrir MATLAB ou Octave.
2. Ajouter le dossier src au chemin :
   addpath(genpath('src'));
3. Lancer un script principal, par exemple :
   principal_temporel

## Références

- Énoncé TP : voir docs/TP_Ondes.pdf
- Rapport : voir docs/Rapport_de_TP_ANN Maria SAAD.pdf

