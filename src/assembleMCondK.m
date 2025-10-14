function [MCond, KK] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembleMCondK :
% assemble les matrices de masse condensée et de raideur globales en P1 lagrange.
%
% SYNOPSIS [MCond, K] = assembleMCondK(Coorneu, Numtri, Reftri)
%          
% INPUT  * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        * Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        * Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%        * Reftri : reference des triangles (matrice entiere Nbtri x 1)
%
% OUTPUT * M matrice de masse globale (vecteur Nbpt)
%        * K matrice de raideur globale (matrice NbptxNbpt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbpt = size(Coorneu, 1);
Nbtri = size(Numtri, 1);

% Declarations des matrices EF.
KK = sparse(Nbpt,Nbpt);
MCondDiag = zeros(Nbpt, 1);

% Boucle d'assemblage sur les triangles.
for l=1:Nbtri
    
    S1 = Coorneu(Numtri(l,1), :);
    S2 = Coorneu(Numtri(l,2), :);
    S3 = Coorneu(Numtri(l,3), :);
    Ref = Reftri(l);
    Kel=matK_elem(S1, S2, S3, Ref);
    % Assemblage de la matrice de rigidité.
    aire_Tl = abs(0.5 * ((S2(1) - S3(1))*(S3(2) - S1(2)) - (S3(1) - S1(1))*(S2(2) - S3(2))));
    
    % Assemblage de la diagonale de la matrice de masse.
    for i=1:3
        MCondDiag(Numtri(l,i))=MCondDiag(Numtri(l,i))+aire_Tl/3;
        for j=1:3
            KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
        end
    end
end % for l

% Transformation de la diagonale en une matrice sparse.
MCond = spdiags(MCondDiag, 0, Nbpt, Nbpt);

% Pseudo élimination.
for l=1:Nbpt
    if or(Refneu(l)==1,Refneu(l)==2)
        MCondDiag(l)=1;
        KK(l,:)=0;
        KK(:,l)=0;
        KK(l,l)=1; %tous les termes sont nuls sauf les diagonaux 
    end
end
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024

