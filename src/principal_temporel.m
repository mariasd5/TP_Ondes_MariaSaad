% =====================================================
% principal_temporel;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante des ondes en regime temporel, avec conditions de
% Dirichlet homogene
% | d^2_{tt} u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

%% Reading mesh file and assembling matrices.

meshFilePath = 'geomRect_0,14.msh';
massType = 'Exacte'; % 'Cond'; %  

[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri]=lecture_msh(meshFilePath);

%Assemblage des matrices M et K 
[M, K] = assembleMCondK(Coorneu, Refneu, Numtri, Reftri);


%% Computing CLF condition and effective time step.

% Calcul de la CFL
cfl=2/(sqrt(abs(max(eigs((M^-1)*K)))));
% Calcul du pas de temps.
cfl_factor = 0.95; %si on prend un delta t plus grand que CFL
dt = cfl_factor * cfl;


%% Computing initial conditions.

interpU0 = zeros(Nbpt, 1);
interpU1 = zeros(Nbpt, 1);
%boucle sur le nombre de sommets
for i = 1:Nbpt
    interpU0(i)=exp(-50*((Coorneu(i,1)-3)^2+(Coorneu(i,2)-1)^2));
end


%% Propagating.

Tmax = 3.0;
niter = floor(Tmax / dt) + 1;
[Us, Kinetic, Potential, Times] = propage_cond(M, K, interpU0, interpU1, dt, niter, Coorneu);



%% Plots of energy and solution.
ouinon=["oui","non"];
plot_energy = "oui";
plot_interp = "oui";
plot_sol = "non";
if (strcmp(plot_energy,"oui"))
    figure
    hold on;
    plot(Times, Potential)
    plot(Times, Kinetic)
    plot(Times, Kinetic + Potential)
    xlim([min(Times) max(Times)])
    xlabel('Time')
    ylabel('Energy')
    legend('P', 'K', 'E', 'Location', 'SouthEast')
    
    %solution a l instant final
    affiche(Us(:, end), Numtri, Coorneu);
    %coefficient sigma
    afficheSigma(Numtri, Reftri, Coorneu);
end

%% Interpolation at point.
if (strcmp(plot_interp,"oui"))
    CoordsInterpPnts = [2, 1; 4, 1];
    NbInterpPnts = size(CoordsInterpPnts, 1);

    % Calcul de la matrice des co?fficients d'interpolation.
    interpolationOp = interpTriP1(Coorneu, Numtri, CoordsInterpPnts);

    s1 = zeros(1, size(Us,2));
s2 = zeros(1, size(Us,2));

for i=1:size(Us,2)
    s1(i) = interpolationOp(1,:)*Us(:,i);
    s2(i) = interpolationOp(2,:)*Us(:,i);
end

end

if (strcmp(plot_sol,"oui"))
for Nt=1:niter
    affiche(Us(:, Nt), Numtri, Coorneu);
    pause(1);
    close all;
end

end
  % Tracé des solutions interpolées aux points (2,1) et (4,1)

figure;
I2_1 = plot([1:length(s1)], s1);
hold on;
I4_1 = plot([1:length(s2)], s2);
legend([I2_1, I4_1], ["Solution en (2,1)", "Solution en (4,1)"]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024

