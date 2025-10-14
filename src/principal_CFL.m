% =====================================================
% principal_CFL;
%
% une routine permettant d'évaluer la dépendance de la condition CFL
% en fonction du pas de maillage
%
% =====================================================

clear 
close all
%% Definition of mesh files and associated mesh steps.

meshFilePath = [ 
    "geomRect_0,02.msh",
    "geomRect_0,08.msh",
    "geomRect_0,14.msh",
    "geomRect_0,2.msh"
    ];

steps = [
    0.02,
    0.08,
    0.14,
    0.2
    ];


%% Computing CFL for every meshes.

nbMesh = length(steps);
if nbMesh ~= length(meshFilePath)
    error('Number of mesh files and number of steps differ.');
end

% Allocation.
cfl = zeros(nbMesh, 1);

% Boucle sur les maillages.
for i=1:nbMesh
    
    [Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(meshFilePath(i));
    
    % Assemblage de M et K et calcul de la CFL.
    [MCond, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri);
    cfl(i)=2/(sqrt(abs(max(eigs((MCond^-1)*K)))));
end


%% Affichage de la cfl en fonction de h.

scatter(steps, cfl);
% Regression lineaire de degre 1
lin = polyfit(steps, cfl, 1);
% Calcule les valeurs 
p = lin(1) .* steps + lin(2) * [1; 1; 1;1];

% bleu pointillé avec des cercles
plot(steps, p, 'b--o');

xlabel('$h$ (Pas de maillage)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\Delta t$ (Condition CFL)', 'Interpreter', 'latex', 'FontSize', 12);

% Annotation pour afficher la pente (coefficient directeur)
coeff_directeur = lin(1); % Récupère la pente
annotation_text = sprintf('Coefficient directeur: %.3f', coeff_directeur);
text(mean(steps), mean(cfl), annotation_text, 'FontSize', 10, 'Color', 'red');

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024


