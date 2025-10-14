function Mel = matM_elem(S1, S2, S3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM_elem :
% calcul la matrices de masse elementaire en P1 lagrange
%
% SYNOPSIS Mel = matM_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture.
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Calcul B_l
B_l = [x2-x1 x3-x1;y2-y1 y3-y1];
% calcul de la matrice de masse.
matrice = ones(3,3) + eye(3,3);
Mel = (abs(det(B_l))/24)*matrice; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024


