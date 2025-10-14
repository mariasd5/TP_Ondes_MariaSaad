function affiche(UU,Numtri,Coorneu, titre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% affiche:
% permet de voir le vecteur UU sur le maillage (Numtri, Coorneu)
%
% SYNOPSIS : affiche(UU,Numtri,Coorneu)
%          
% INPUT * UU vecteur de valeurs aux sommets (vecteur Nbpt x 1)
%       * Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 3)
%       * Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%       * titre (optionel) un titre (string)
%
% OUTPUT une fenetre graphique
%
% NOTE 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy = max(Coorneu(:,2)) - min(Coorneu(:,2));
dx = max(Coorneu(:,1)) - min(Coorneu(:,1));
ratio = dx / dy;

figure; 
axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
            -1,1,-1,1]);
% control on the input args
if (nargin<4) 
    titre = ''; 
end
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),UU);
pbaspect([ratio 1 1])
view(2);
shading interp
% shading faceted
% shading flat
colorbar;

title("visualisation du vecteur UU sur le maillage");
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024



