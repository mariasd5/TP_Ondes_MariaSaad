function Kel = matK_elem(S1, S2, S3, Reftri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matK_elem :
% calcul la matrices de raideur elementaire en P1 lagrange.
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%       * Reftri : reference du triangle.
%
% OUTPUT * Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul utilise une formule de quadrature de Gauss-Legendre.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Reftri ~= 1 && Reftri ~= 2
    error('Un des triangles a une réference differente de 1 et 2.');
end

% preliminaires, pour faciliter la lecture.
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Points et poids de quadrature.
S_hat(:,1) = [1/6; 1/6];
S_hat(:,2) = [2/3; 1/6];
S_hat(:,3) = [1/6; 2/3];
w0 = 1/6;

% Gradients des fonctions de base sur le triangle de reference.
grad_w1 = [-1;-1];
grad_w2 = [1;0];
grad_w3 = [0;1];
grad_w = [grad_w1 grad_w2 grad_w3];


% Transformation géométrique associée au triangle courant.
 B_l = [x2-x1 x3-x1;y2-y1 y3-y1];

% Boucle sur les fonctions de bases locales.
% Calcul de l'integrale de sigma a l'aide de Gauss-Legendre
% Il faut verifier si c'est dans la zone 1 ou 2 a l'aide de reftri
if Reftri == 1
    int_sigma = w0*(sigma_1(S_hat(1,1),S_hat(2,1))+sigma_1(S_hat(1,2),S_hat(2,2))+sigma_1(S_hat(1,3),S_hat(2,3)));
else  %ie Reftri=2
    int_sigma = w0*(sigma_2(S_hat(1,1),S_hat(2,1))+sigma_2(S_hat(1,2),S_hat(2,2))+sigma_2(S_hat(1,3),S_hat(2,3)));
end

%Calcul matrice de raideur elementaire
Kel = zeros(3,3);
for j=1:3
    for i=1:3
        Kel(i,j)=abs(det(B_l))*int_sigma*(((((B_l)')^-1)*grad_w(:,i))')*((((B_l)')^-1)*grad_w(:,j));
    end
end



