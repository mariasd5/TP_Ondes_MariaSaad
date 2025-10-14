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
Us(:, 1) = % A COMPLETER
Us(:, 2) = % A COMPLETER

% Iteriations.
for i = 1:niter
    
    % Calcul des ?nergies.
    % A COMPLETER
    
    % Calcul de la solution ? l'it?ration suivante par descente remont?e.
    % A COMPLETER
    
    Times(i) = i * dt;
end

end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024

