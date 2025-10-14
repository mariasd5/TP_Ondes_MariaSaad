function val = f(x, y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val =(1+(1/81+1/4)*(pi^2))*sin(pi*x/9)*sin(pi*y/2); %1.12
val =600*exp(-((2*x-9)/7.2)^2-((y-1)/0.8)^2) ;%1.15

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2024


