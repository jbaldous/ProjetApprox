//Modèle à 1 ressort : 

//Paramètres :

E = 200*10^9 //Module de Young
m = 800 //Masse de la pale
rho = 8000 //Densité de la pale
L = 14 // Longeur de la pale
S = m/(rho*L) //Section équilvalente de la pale
Omega0 = 2*%pi //Vitesse angulaire de la pale
k = E*S/L //Raideur du ressort

xeq = L + m*(Omega0^2)*L/(2*k) //Solution stationnaire du modèle

xini = L //Condition initiale pour x(0)
vini = L*Omega0 //Vitesse initiale pour x'(0)

//Autres conditions initiales intéressantes :

//xini = xeq
//vini = 0

//xini = xeq + (rand(1)-rand(1)) * L/100
//vini = 0

T = 0.5 //Temps max de la subdivision
N = 100000 //Nombre de points de la subdivision
t = linspace(0,T,N+1) //Intervalle de la subdivision

//Fonctions :

//Fonction pour calculer les approximations de la solution exacte par Euler Explicite
function X = Euler_exp(T,N,k,m,L,Omega0,xini,vini)
    delta = T/N //Pas de la subdivision
    X = zeros(2,N+1)
    X(:,1) = [xini;vini] //Conditions initiales du schéma numérique
    
    for n = 1 : N //Boucle calculant les itérées du schéma numérique
        X(1,n+1) = X(1,n) + delta * X(2,n)
        X(2,n+1) = X(2,n) + delta * (- 2*k/m * X(1,n) + 2*k*L/m + (Omega0^2) * L)
    end
endfunction

//Fonction pour calculer les approximations de la solution exacte par Cranck-Nicholson
function X = Cranck(T,N,k,m,L,Omega0,xini,vini)
    delta = T/N //Pas de la subdivision
    X = zeros(2,N+1)
    X(:,1) = [xini;vini] //Conditions initiales du schéma numérique
    A = [0,1;-2*k/m,0] //Matrice pour simplifier les notations
    
    for n = 1 : N //Boucle calculant les itérées du schéma numérique
        X(:,n+1) = inv(diag(ones(1,2)) - delta/2*A) * ((diag(ones(1,2)) + delta/2*A) * X(:,n) + delta * [0;2*k*L/m + Omega0^2*L] )
    end
endfunction


//Fonction pour calculer la solution exacte du modèle
function x = sol_exacte(T,N,k,m,L,Omega0,xini,vini,xeq)
    x = xeq + (xini-xeq)*cos(sqrt(2*k/m)*t) + sqrt(m/(2*k))*vini*sin(sqrt(2*k/m)*t)
endfunction


//Représentation graphique :

//Figure 1
scf(0)

//Représentation des vecteurs calculés
plot2d(t',[sol_exacte(T,N,k,m,L,Omega0,xini,vini,xeq);Euler_exp(T,N,k,m,L,Omega0,xini,vini)(1,:)]',[2, 33]);

//Légende
legends(['Solution Exacte';'Euler explicite'],[2,33],opt="lr")

//Titre
title('$Résolution \hspace{0.1cm} numérique \hspace{0.1cm} du  \hspace{0.1cm} P. C.  \hspace{0.1cm} avec \hspace{0.1cm} x(0) = L \hspace{0.1cm} et \hspace{0.1cm} x''(0) = \Omega_0L$')

//Grille
gcf().children.grid = color("grey70")*[1 1]


//Figure 2
scf(1)

//Représentation des vecteurs calculés
plot2d(t',[sol_exacte(T,N,k,m,L,Omega0,xini,vini,xeq);Cranck(T,N,k,m,L,Omega0,xini,vini)(1,:)]',[2, 33]);

//Légende
legends(['Solution Exacte';'Cranck-Nicholson'],[2,33],opt="lr")

//Titre
title('$Résolution \hspace{0.1cm} numérique \hspace{0.1cm} du  \hspace{0.1cm} P. C.  \hspace{0.1cm} avec \hspace{0.1cm} x(0) = L \hspace{0.1cm} et \hspace{0.1cm} x''(0) = \Omega_0L$')

//Grille
gcf().children.grid = color("grey70")*[1 1]
