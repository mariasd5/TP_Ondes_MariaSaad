% =====================================================
% principal_stationnaire;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante stationnaire, avec conditions de
% Dirichlet homogene
% | u - div(\sigma \grad u)= f,   dans \Omega=\Omega_1 U \Omega_2
% |         u = 0,   sur le bord
%
% avec
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% =====================================================

% Donnees du probleme.
nom_maillage = 'geomRect.msh';
affichage = true; % false; %

% Lecture du maillage et affichage.
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri] = lecture_msh(nom_maillage);

% Declarations des matrices EF et vecteur second membre.
KK = sparse(Nbpt,Nbpt);  %car les matrices sont creuses 
MM = sparse(Nbpt,Nbpt);
FF = zeros(Nbpt,1);

% Boucle d'assemblage.
% Boucle sur les triangles 
for l=1:Nbtri
    % On recupere les sommets du triangle
    S1 = Coorneu(Numtri(l,1),:);
    S2 = Coorneu(Numtri(l,2),:);
    S3 = Coorneu(Numtri(l,3),:);
    % Pour savoir si on est dans omega 1 ou 2
    Ref = Reftri(l);
    % Calcul des matrices elementaires du triangle l.
    Kel = matK_elem(S1, S2, S3, Ref);
    Mel = matM_elem(S1, S2, S3);
    
    % Assemblage des matrices globales.
  for i=1:3
        for j=1:3
            MM(Numtri(l,i),Numtri(l,j))=MM(Numtri(l,i),Numtri(l,j))+Mel(i,j);
            KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
        end
  end    
end

% Calcul du second membre.
for i = 1:Nbpt
    FF(i) = f(Coorneu(i,1),Coorneu(i,2));
end
LL = MM*FF;
% Matrice EF complete.
AA = MM + KK;

% Pseudo-elimination.
[tilde_AA, tilde_LL] = elimine(AA,LL,Refneu);

% R?solution du prolb?me par inversion.
UU = tilde_AA\tilde_LL;
maxValue = max(UU)
% Visualisation de sigma et de la solution.
if affichage
    afficheSigma(Numtri, Reftri, Coorneu);
    affiche(UU, Numtri, Coorneu, sprintf('Stationnaire - %s', nom_maillage));
end

 %tracer de l'erreur:
 figure;
hold on
axis('equal');
title('Régression linéaire pour h=0.5 0.2 0.05');
h=[log(0.5); log(0.2); log(0.05)];

% Calcul de la matruice U exact
% Récupération des coordonnées x et y
x = Coorneu(:, 1);
y = Coorneu(:, 2);
U_ex = sin(pi / 9 * x) .* sin(pi / 2 * y);
erreur_L2 = [-1.59976 -2.3371 -3.5338]; %calcule a la main dans un fichier excel 
erreur_H1 = [-1.4859 -2.1934 -3.2778];  %calcule a la main dans un fichier excel 
%NUM_L2 = sqrt((U_ex-UU)'*MM*(U_ex-UU))
%NUM_H1 = sqrt((U_ex-UU)'*KK*(U_ex-UU))
%DENOM_L2=sqrt((U_ex)'*MM*(U_ex))
%DENOM_H1=sqrt((U_ex)'*KK*(U_ex))

plot(h,erreur_H1,'r');
plot(h,erreur_L2,'b');

pL2 = polyfit(h,erreur_L2,1);
pH1 = polyfit(h,erreur_H1,1);

p1=pH1(1).*h + pH1(2)*[1;1;1];
p2=pL2(1).*h + pL2(2)*[1;1;1];

plot(h,p1,'r--o');
plot(h,p2,'b--o');
xlabel('log(h)');       % Nom de l'axe des abscisses
ylabel('log(erreur)');  % Nom de l'axe des ordonnées

% Ajout d'une légende
legend('H1 (rouge)', 'L2 (bleu)', 'Location', 'best');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024



