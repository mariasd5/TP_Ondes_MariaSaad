function [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter, Coorneu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage_cond :
% Propage les conditions initiales à partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter)
%
% INPUT  * M, K : matrice de masse condensée et de rigidité
%        * interpU0, interpU1 : interpolée en espace des conditions initiales
%        * dt : pas de temps du schéma.
%        * niter : nombre d'itérations
%
% OUTPUT * Us : saolution discrète à tous les pas de temps
%               avec les conditions initiales (matrice Nbpt x niter + 2)
%        * Kinetic, Potential : énergies cinétique et potentielle (vecteur niter)
%        * Times : vecteur temps (vecteur niter)
%
% NOTE (1) On suppose que la matrice de masse est diagonale de sorte que
% l'opération M \ b est efficace.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation.
Nbpt = length(interpU0);
Us = zeros(Nbpt, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);

% Conditions initiales.
for i =1:Nbpt
    Us(i, 1) = exp(-50*((Coorneu(i,1)-3)^2+(Coorneu(i,2)-1)^2));
end
Us(:, 2) = zeros(Nbpt,1);

% Interiations.
for i = 1:niter
    
    % Calcul des énergies.
    Kinetic(i) = 0.5 * ((Us(:, i+1) - Us(:, i)) / dt)' * (M - (dt^2 / 4) * K) * ((Us(:, i+1) - Us(:, i)) / dt);
    Potential(i) = 0.5 * ((Us(:, i+1) + Us(:, i)) / 2)' * K * ((Us(:, i+1) + Us(:, i)) / 2);

    
    % Calcul de la solution par résolution directe du système linéaire
    % (suppose que l'opération M \ b soit efficace, i.e. M diagonale).
    Us(:, i+2) = M \ ((2 * M - (dt^2) * K) * Us(:, i+1) - M * Us(:, i));
    Times(i) = i * dt;
end
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024

