function [Us, Kinetic, Potential, Times] = propage_implicite(M, K, interpU0, interpU1, dt, niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propage :
% Propage les conditions initiales ? partir du propagateur discret saute-mouton.
%
% SYNOPSIS [Us, Kinetic, Potential, Times] = propage(M, K, interpU0, interpU1, dt, niter)
%
% INPUT  * M, K : matrice de masse et de rigidit?
%        * interpU0, interpU1 : interpol?e en espace des conditions initiales
%        * dt : pas de temps du sch?ma.
%        * niter : nombre d'it?rations
%
% OUTPUT * Us : saolution discr?te ? tous les pas de temps
%               avec les conditions initiales (matrice Nbpt x niter + 2)
%        * Kinetic, Potential : ?nergies cin?tique et potentielle (vecteur niter)
%        * Times : vecteur temps (vecteur niter)
%
% NOTE (1) On utilise la d?composition de Cholesky pour r?soudre les
% syt?mes du type M X = b dans les it?rations en temps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul de la d?composition de Cholesky de la matrice de masse.
cholU = chol(M+dt^2/4*K);
% pour calculer (M+dt^2/4*K)^{-1}*V  il suffit d ecrire cholU \ (cholU' \ V);

% Allocation.
Nbpts = length(interpU0);
Us = zeros(Nbpts, niter + 2);
Kinetic = zeros(niter, 1);
Potential = zeros(niter, 1);
Times = zeros(niter, 1);

% Conditions initiales.
Us(:, 1) = interpU0;
Us(:, 2) = interpU0 + dt * interpU1 - (dt^2 / 2) * (M \ (K * interpU0));


% Iteriations.
for i = 1:niter
    
    % Calcul des énergies
    vitesse = (Us(:, i) - Us(:, i-1)) / dt;
    deplacement = (Us(:, i) + Us(:, i-1)) / 2;
    Kinetic(i) = 0.5 * (vitesse' * (M * vitesse));
    Potential(i) = 0.5 * (deplacement' * (K * deplacement));

    % Calcul de la solution suivante (schéma implicite)
    rhs = (2 * M - (dt^2 / 2) * K) * Us(:, i) - (M - (dt^2 / 4) * K) * Us(:, i-1);
    Us(:, i+1) = (M + (dt^2 / 4) * K) \ rhs;

    
    Times(i) = i * dt;
end

end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024

