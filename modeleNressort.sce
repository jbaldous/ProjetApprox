//Modèle à N ressorts : 
//Paramètres :

N = 10 //Nombres de ressorts

E = 200*10^9 //Module de Young
m = 800 //Masse de la pale
rho = 8000 //Densité de la pale
L = 14 // Longeur de la pale
S = m/(rho*L) //Section équilvalente de la pale
Omega0 = 2*%pi //Vitesse angulaire de la pale
k = E*S/L //Raideur du ressort

T = 0.5 //Temps max de la subdivision
n = 100000 //Nombre de points de la subdivision
t = linspace(0,T,n+1) //Intervalle de la subdivision

//Matrices et vecteurs : 
//Matrice du P. C. :
A1 = -2 * eye(N,N) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1)
A1(N,N) = -1
A1(N,N-1) = 1

//Vecteur du P. C. :
v1 = zeros(N,1) //Vecteur du P. C.
v1(N) = m * Omega0^2 * L/(2*N) + k*L/N
for i = 1 : N-1
    v1(i) = m * i * (Omega0/N)^2 * L
end

//Solution stationnaire du modèle :
//Xeq = -inv(k*A1) * v1

//Matrice A pour le schéma :
A1_mod = k*N/m*A1
A1_mod(N,N) = 2*A1_mod(N,N)
A1_mod(N,N-1) = 2*A1_mod(N,N-1)

A = [zeros(N,N),  diag(ones(N,1)) ; A1_mod , zeros(N,N)] //Matrice du schéma de Cranck-Nicholson

//Vecteur F pour le schéma :
F = zeros(2*N,1)
for i = N+1 : 2*N-1
    F(i) =  i * Omega0^2 * L/N
end
F(2*N) = Omega0^2 * L + 2*k*L/m


//Solution à l'équilibre du problème à l'orde 1 :
Yeq = -inv(A)*F

//On en déduit la solution à l'équilibre du P. C. (ordre 2) :
Xeq = Yeq(1:N)

//Conditions initiales :

//Xini = Xeq
//Vini = Yeq(N+1:2*N)

//Autres conditions initiales intéressantes :

Xini = Xeq + (rand(1)-rand(1)) * (Xeq(N) - L)
Vini = Yeq(N+1:2*N)

//Fonctions :

//Fonction pour calculer les approximations de la solution exacte par Cranck-Nicholson
function X = Cranck(N,T,n,Xini,Vini)
    delta = T/n //Pas de la subdivision
    X = zeros(2*N,n+1)
    X(:,1) = [Xini;Vini] //Conditions initiales du schéma numérique
    
    for i = 1 : n //Boucle calculant les itérées du schéma numérique
        X(:,i+1) = inv(diag(ones(1,2*N)) - delta/2*A) * ((diag(ones(1,2 * N)) + delta/2*A) * X(:,i) + delta * F )
    end
endfunction

//Calcul de toutes les positions des masses au temps final (en incluant x0) :
Xn = [0;Cranck(N,T,n,Xini,Vini)(1:N,n+1)]

//Calcul de l'élongation par rapport à la longueur de la pale L :
Xn = Xn - [0:L/N:L]'

//Représentation graphique :
//Figure 0
scf(0)

//Représentation des vecteurs calculés
plot2d([0:N],Xn);

//Titre
title('$Elongation \hspace{0.1cm} des \hspace{0.1cm} masses  \hspace{0.1cm} avec  \hspace{0.1cm} 10 \hspace{0.2cm} ressorts \hspace{0.1cm} au \hspace{0.1cm} temps \hspace{0.1cm} final$')

//Grille
gcf().children.grid = color("grey70")*[1 1]

//Figure 1
scf(1)

//Représentation des vecteurs calculés
plot2d(t',Cranck(N,T,n,Xini,Vini)(N,:)',3);
plot2d(t',Xeq(N)*ones(n+1,1),7);

//Légende
legends(['Solution stationnaire';'Cranck-Nicholson avec condition à l''équilibre perturbé'],[7,3],opt="?")

//Titre
title('$Résolution \hspace{0.1cm} numérique \hspace{0.1cm} du  \hspace{0.1cm} P. C.  \hspace{0.1cm} avec \hspace{0.2cm} 10 \hspace{0.2cm} ressorts \hspace{0.1cm}$')

//Grille
gcf().children.grid = color("grey70")*[1 1]
