
% ANALYSE DE LA RUPTURE D'UNE PLAQUE
% COMPOSITE EN TRACTION/FLEXION
% TP MATLAB DESS UFR MIG
% Version 2018

% 07 decembre 2018
% POUR TP DE 3h
% Ajouter les case
% corriger le flambement 
% Ajouter les contraintes cas MEMBRANE
% Ajouter le cas epaisseur avec pourcentage
% Ajouter la flexion et calcul des contraintes..
%
% Pour le projet faire le Sandwich en flexion
% Choix flambement encastree

% Initialisation

clear all;
close all;

% UNITES 
% mm MPa s kg

% Données matériaux commune au pli
El=165000;      % module sens fibres
Et=8500;        % module sens travers E22
Etp=8900;       % module sens travers hors plan E33
Glt=4800;
Gltp=5300;      % G13
Gttp=6700;      % G23
nult=0.31;      % 
nultp=0.31;
nuttp=0.35;
nutl=nult*Et/El;
alpha=1/(1-nult*nutl);
ep=0.25;            % en mm

% Valeur à rupture (MPa)
siglrt=2600;siglrc=1500;sigtrt=80;sigtrc=400;sigltr=95;
% Dimension de la plaque rectangulaire
% a = longueur , b = largeur
a=200;b=10;k=0;


disp('Choix 1 : Raideur globale d un pli en fonction de theta')
disp('Choix 2 : Effort maximal admissible d un stratifie en membrane')

disp('Choix 3 : Epaisseur minimale d un stratifie en membrane')
disp('Choix 4 : Dimensionnement d un stratifié en flexion')
disp('Choix 5 : Plaque orthotrope appuyee en flambement bi-axiale')
disp('Choix 6 : calcul d une poutre sandwich en flexion')
disp('other value')

n = input('Entrez votre choix : ');

switch n

%==========================================================================
    case 1
        disp('raideur d un pli seul en fonction de angle')
        
% Matrice de rigidité du pli dans le repère local
Slt=[1/El -nutl/Et 0 ; -nult/El 1/Et 0 ; 0 0 1/Glt];
Qlt=inv(Slt);

% Declaration des tableaux
Qxy=zeros(3);
Sxy=zeros(3);
epsx=zeros(3,1);
% Pour critere de Tsai-Hill
sigl=zeros(3,91);
crit1=zeros(1,91);
sigxm1=zeros(1,91);
% Chargement
sigx=[1,0,0];
sigx=sigx';

% Si critere de Hashin
crit2=zeros(4,91);
sigxm2=zeros(4,91);

% % Valeur à rupture (MPa)
% siglrt=2600;siglrc=1400;sigtrt=80;sigtrc=400;sigltr=95;


for i=0:90
    c=cos(i*pi/180);
    s=sin(i*pi/180);
     
    TS=[c*c s*s  -2*s*c ; s*s c*c  2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s  -s*c ; s*s c*c  s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'
     
% Matrice de rigidité/souplesse du pli dans le repère global

    Qxy=TS*Qlt*inv(TE);
    Sxy=inv(Qxy);                                            % ou Sxy=TE*Slt*inv(TS);
    
    epsx=Sxy*sigx;
    sigl(:,i+1)=Qlt*inv(TE)*epsx;
    
    if sigl(1,i+1)>=0
        srl=siglrt;
    else
        srl=siglrc;
    end
    if sigl(2,i+1)>=0
        srt=sigtrt;
    else
        srt=sigtrc;
    end
        crit1(1,i+1)=(sigl(1,i+1)/srl)^2+(sigl(2,i+1)/srt)^2-(sigl(1,i+1)*sigl(2,i+1)/(srl^2))+(sigl(3,i+1)/sigltr)^2;
        sigxm1(1,i+1)=1/sqrt(crit1(1,i+1));  
        
        
% Module et coefficient de Poisson dans le repère global
    Ex=1/Sxy(1,1);
    Ey=1/Sxy(2,2);
    Gxy=1/Sxy(3,3);
    
    nuxy=-Ex*Sxy(2,1);
    nuyx=-Ey*Sxy(1,2);
        
    muxy=-Ey*Sxy(3,2);
    muyx=-Gxy*Sxy(2,3);
    
    etaxy=-Ex*Sxy(3,1);
    etayx=-Gxy*Sxy(1,3);

% Sauvegarde dans des tableaux pour graphique
Exx(i+1,1)=Ex;
Eyy(i+1,1)=Ey;
Gxxy(i+1,1)=10.*Gxy;           % Gxy *10 pour le tracé

nuxxy(i+1,1)=nuxy;
nuxyx(i+1,1)=nuyx;
etaxxy(i+1,1)=etaxy;
etaxyx(i+1,1)=etayx;
muxxy(i+1,1)=muxy;
muxyx(i+1,1)=muyx;
end

% tracé des modules et des coef de Poisson
%==========================================
figure
subplot(2,2,1)
plot(0:90,Exx,0:90,Eyy,0:90,Gxxy)
xlabel('theta (°)')
ylabel('Module (MPa)')
%axis([0 90 0 150000])
grid on
subplot(2,2,2)
plot(0:90,nuxxy,0:90,nuxyx)
%axis([0 90 0 0.5])
grid on
subplot(2,2,3)
plot(0:90,muxxy,0:90,muxyx)
%axis([0 90 0 3])
grid on
subplot(2,2,4)
plot(0:90,etaxxy,0:90,etaxyx)
axis([0 90 0 3])
grid on

% Tracé de la contrainte max à rupture
%======================================
figure
plot(0:90,sigxm1)
grid on
    
    
    case 2
 disp('Effort maximal admissible d un stratifie en membrane')
%========================================================================
% Calcul à la rupture d'une séquence d'empilement
%========================================================================
% init = sequence d'empilement
% CALCUL DE LA MATRICE de COMPORTEMNT COMPLETE 
% Calcul des modules homogénéisés de la séquence

% données matériaux séquences....
%seq=[0, -45 ,90,90,-45 ,0]; % en degré
%seq=[45, 0 ,-45,90,90 ,-45,0,45];
%seq=[45,0,-45,90,45,0,-45,90,90,90,90,-45,0,45,90,-45,0,45]; % en degré
%seq=[0,45,90,-45,0,45,90,-45,-45,-45,-45,90,45,0,-45,90,45,0]; % en degré
%seq=[0,90,90,0];


% Panneau SUP et INF
%seq=[45,0,0,-45,0,0,90,90,0,0,-45,0,0,45];

%seq=[45,-45,0,0,90,0,0,90,0,0,90,0,0,-45,45]

seq=[0,0,45,0,0,-45,90,0,90,-45,0,0,45,0,0];

% Panneau Lateraaux
%seq=[45, -45 ,0,90,0 ,-45,45];
%seq=[45, -45 ,0,90,90,0 ,-45,45];



nbp=length(seq);

ht=ep*length(seq);

% si calcul de epaisseur totale
p=[0.25,0.25,0.25,0.25];
hc=p(1)+p(2)+p(3)+p(4);

% Chargement (N/mm)
disp('chargement de type Nxx = 1 ou 0 ; Nyy = 1 ou 0 ; Nxy = 1 ou 0')
Nx=zeros(3,1);
Nx = input('Entrer le chargement en membrane : '); 
%Nx=[1,0,0];
Nx=Nx';
% % Valeur à rupture (MPa)
% siglrt=2600;siglrc=1500;sigtrt=80;sigtrc=400;sigltr=95;

% initialisation des tableaux
A=zeros(3);
B=zeros(3);
D=zeros(3);
epsx=zeros(3,1);
sigl=zeros(3,nbp);
crit=zeros(1,nbp);
Nmax=zeros(1,nbp);

% Matrice de rigidité du pli dans le repère local
Slt=[1/El -nutl/Et 0 ; -nult/El 1/Et 0 ; 0 0 1/Glt];
Qlt=inv(Slt);

for i=1:(nbp+1)
    h(i)=-ht/2+(i-1)*ep;  
end

for i=1:nbp
    c=cos(seq(i)*pi/180);
    s=sin(seq(i)*pi/180);
    TS=[c*c s*s -2*s*c ; s*s c*c 2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s -s*c ; s*s c*c s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'     
% Matrice de rigidité/souplesse du pli dans le repère global

    Qxy=TS*Qlt*inv(TE);
    
    A=A+Qxy*(h(i+1)-h(i));
    B=B+(1/2)*Qxy*(h(i+1)^2-h(i)^2);
    D=D+(1/3)*Qxy*(h(i+1)^3-h(i)^3);
end
% Module Homogénéisé
   A1=inv(A);
   D1=inv(D);
ABD=[A B; B D];
ABD_I=inv(ABD);

%ABD



% Module de traction  
disp('module global Ex du stratifié');
   Exg=1/(ht*A1(1,1));
   Eyg=1/(ht*A1(2,2));
   Ezg=1/(ht*A1(3,3));
   Exgd=1/(ht*ABD_I(1,1));
   
   Gxyg=1/(ht*A1(3,3));
%   Gxzg=1/(ht*A1(5,5));
   
   nuxyg=-ht*Exg*A1(1,2);
   nuyxg=-ht*Eyg*A1(2,1);
   Exg,Eyg,Gxyg,nuxyg %,Gxyg,Gxzg,nuxyg,nuyxg
   
% Module apparent de Flexion 
   Efx=1/((ht^3/12)*D1(1,1));
   Efy=1/((ht^3/12)*D1(2,2)); 
   Efx,Efy
    
% Déformation globale du stratifié en traction 
   epsx=A1*Nx; 
   
% calcul des contraintes dans les couches
for i=1:nbp;
    c=cos(seq(i)*pi/180);
    s=sin(seq(i)*pi/180);
    TS=[c*c s*s -2*s*c ; s*s c*c 2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s -s*c ; s*s c*c s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'
    
    sigl(:,i)=Qlt*inv(TE)*epsx;
    
    if sigl(1,i)>=0
        srl=siglrt;
    else
        srl=siglrc;
    end
    if sigl(2,i)>=0
        srt=sigtrt;
    else
        srt=sigtrc;
    end
        crit(1,i)=(sigl(1,i)/srl)^2+(sigl(2,i)/srt)^2-(sigl(1,i)*sigl(2,i)/(srl^2))+(sigl(3,i)/sigltr)^2;
        Nmax(1,i)=1/sqrt(crit(1,i));
% calcul du critere de Tsai-Hill
end
    
    sigl    
    Nmax
    
% Calcul en flambement globale
% =========================================================================
% Version Sandwich H
    H=18.;b=50;
    
% Effort critique de flambement initial  
    disp('Nx critique initial simplement appuy�e');
    Ncr_init=2.*3.14115^2/b^2*(sqrt(D(1,1)*D(2,2))+(D(1,2)+2.*D(3,3)))
    
    disp('Nx critique initial encastr�e');
    Ncr_init=1./b^2*(44.6*sqrt(D(1,1)*D(2,2))+2.46*3.14115^2*(D(1,2)+2.*D(3,3)))
    
% Effort critique de flambement initial  
    D_F=2.*A*((H/2-ht/2)^3-(H/2)^3)/3;
   
    disp('Nx critique sandwich simplement appuy�e');
    Ncr=2.*3.14115^2/b^2*(sqrt(D_F(1,1)*D_F(2,2))+(D_F(1,2)+2.*D_F(3,3)))
    
    disp('Nx critique sandwich encastr�e');
    Ncr_init=1./b^2*(44.6*sqrt(D_F(1,1)*D_F(2,2))+2.46*3.14115^2*(D_F(1,2)+2.*D_F(3,3)))
    
    
%     Nmini=min(Nmax)
%     Nmaxi=max(Nmax)
%Nmn=zeros(1,2500);xx=zeros(1,2500);

%tracé de sigma sens fibres
figure(1);
YC=[length(seq):1:1];
XC=[1];
clim=[min(sigl),max(sigl)];
plot(1,length(seq));
hold on;
image(XC,YC,sigl(1,:)','CDataMapping','scaled');
title('contrainte sens fibres');
colormap(jet);
colorbar

%tracé de sigma sens travers
figure(2);
YC=[length(seq):1:1];
XC=[1 ];
clim=[min(sigl),max(sigl)];
plot(1,length(seq));
hold on;
image(XC,YC,sigl(2,:)','CDataMapping','scaled');
title('contrainte sens travers');
ylabel('sequence d empilement');
colormap(jet);
colorbar

%tracé de sigma cisaillement plan 12
figure(3);
YC=[length(seq):1:1];
XC=[1 ];
clim=[min(sigl),max(sigl)];
plot(1,length(seq));
hold on;

image(XC,YC,sigl(3,:)','CDataMapping','scaled');
title('contrainte de cisaillement');
%im.AlphaData = 0.1;
colormap(jet);
colorbar

    case 3
        
        disp('Epaisseur minimale d un stratifie en membrane')
%========================================================================
% Calcul à la rupture d'une séquence d'empilement
%========================================================================
% init = sequence d'empilement
% CALCUL DE LA MATRICE de COMPORTEMNT COMPLETE 
% Calcul des modules homogénéisés de la séquence

% Chargement (N/mm)

disp('Chargement sous la forme [Nx Ny Nxy]')
Nx=zeros(3,1);
Nx = input('Entrer chargement : '); 

disp('Pourcentage 0.25 = 25% sous la forme [%0 %45 %-45 %90]')
p=zeros(4,1);
theta=[0,45,-45,90];
p = input('Entrer le pourcentage : '); 
Nx=Nx';
hc=p(1)+p(2)+p(3)+p(4);
if (hc>1) 
    disp('Attention pourcentage supérieur à 1')
end
nbp=length(p);

% initialisation des tableaux
A=zeros(3);
B=zeros(3);
D=zeros(3);
epsx=zeros(3,1);
sigl=zeros(3,nbp);
crit=zeros(1,nbp);
emin=zeros(1,nbp);

% Matrice de rigidité du pli dans le repère local
Slt=[1/El -nutl/Et 0 ; -nult/El 1/Et 0 ; 0 0 1/Glt];
Qlt=inv(Slt);

for i=1:nbp
    c=cos(theta(i)*pi/180);
    s=sin(theta(i)*pi/180);
    TS=[c*c s*s -2*s*c ; s*s c*c 2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s -s*c ; s*s c*c s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'     
% Matrice de rigidité/souplesse du pli dans le repère global

    Qxy=TS*Qlt*inv(TE);
    
    A=A+Qxy*p(i);     % A=A+Qxy*ek/ht  ou A=A+Qxy*(hk-hk-1)
    
end
% Module Homogénéisé
   A1=inv(A);
   
% Déformation globale du stratifié en traction 
   epsx=A1*Nx; 
   
% calcul des contraintes dans les couches
for i=1:nbp
    c=cos(theta(i)*pi/180);
    s=sin(theta(i)*pi/180);
    TS=[c*c s*s -2*s*c ; s*s c*c 2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s -s*c ; s*s c*c s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'
    
    sigl(:,i)=Qlt*inv(TE)*epsx;
    
    if sigl(1,i)>=0
        srl=siglrt;
    else
        srl=siglrc;
    end
    if sigl(2,i)>=0
        srt=sigtrt;
    else
        srt=sigtrc;
    end
    % critere de Tsai-Hill
        crit(1,i)=(sigl(1,i)/srl)^2+(sigl(2,i)/srt)^2-(sigl(1,i)*sigl(2,i)/(srl^2))+(sigl(3,i)/sigltr)^2;
        
        emin(1,i)=sqrt(crit(1,i));
        
% calcul du critere de Tsai-Hill

end

    disp(' ');
    disp('===================================================================');
    disp('Contrainte Sig_l, Sig_t, Sig_lt dans les plis � 0�, 45�, -45�, 90� ');
    sigl   

    
    disp('===================================================================');
    disp('Epaisseur minimale pour les plis � 0�, 45�, -45�, 90� ');
    emin
    
%     Nmini=min(Nmax)
%     Nmaxi=max(Nmax)
%Nmn=zeros(1,2500);xx=zeros(1,2500);
        
    case 4
    
     disp ('D�form�e et contraintes dans une Plaque appuy�e en flexion')   

seq=[45,0,90,-45,-45,90,0,45];      % Sequence d'empilement
ep=0.25;                            % Epaisseur du pli en mm
F=-1000;                            % force appliquer en N
effort=2;                           % 1 Effort ponctuel ; 2 pression sur surface c x d
a=400;                              % longeur de la plaque en mm
b=200;                              % largeur de la plaque en mm
c=4;d=2;                            % longueur et largueur du rectangle c/x d/y ou la pression est appliqué
x0=a/4;y0=b/2;                      % Position du centre du rectangle où la force est appliquée

if (((x0-c/2)<0) || ((x0+c/2)>a) || ((y0-d/2)<0) || ((y0+d/2)>b))
   fprintf('\n Attention : Mauvais positionnement de l''effort');
   return;
end
q0=F/(c*d);
nm_max=13;    % Paramètre m et n de la décomposition en série de fourrier


nbp=length(seq);
ht=ep*length(seq);
k=0;
% initialisation des tableaux
%==========================================================================
% CALCUL STRATIFIE
%==========================================================================
A=zeros(3);
B=zeros(3);
D=zeros(3);
epsx=zeros(3,1);
sigl=zeros(3,nbp);
crit=zeros(1,nbp);
Nmax=zeros(1,nbp);

% Matrice de rigidité du pli dans le repère local
Slt=[1/El -nutl/Et 0 ; -nult/El 1/Et 0 ; 0 0 1/Glt];
Qlt=inv(Slt);

for i=1:(nbp+1)
    h(i)=-ht/2+(i-1)*ep;  
end

for i=1:nbp
    c=cos(seq(i)*pi/180);
    s=sin(seq(i)*pi/180);
    TS=[c*c s*s -2*s*c ; s*s c*c 2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s -s*c ; s*s c*c s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'     
% Matrice de rigidité/souplesse du pli dans le repère global

    Qxy=TS*Qlt*inv(TE);
    
    A=A+Qxy*(h(i+1)-h(i));
    B=B+(1/2)*Qxy*(h(i+1)^2-h(i)^2);
    D=D+(1/3)*Qxy*(h(i+1)^3-h(i)^3);
end

ABD=[A B; B D];
% Affichage matrice de raideur du stratifié
ABD
nbpx=50;nbpy=50;      % nombre de points de calculs de w / x et y
w=zeros(nbpx,nbpy);
x=0;
for i=1:nbpx
    x=x+a/nbpx;
    y=0;
    for j=1:nbpy
       y=y+b/nbpy;      
        for m=1:1:nm_max
            for n=1:1:nm_max
    
    Dmn=D(1,1)*(m/a)^4+2*(D(1,2)+2*D(3,3))*((m/a)^2*(n/b)^2)+D(2,2)*(n/b)^4;
    
            if (effort == 1) 
    w(i,j)=w(i,j)+(4*F/(pi^4*Dmn*a*b))*(sin(m*pi*x0/a)*sin(n*pi*y0/b))*sin(m*pi*x/a)*sin(n*pi*y/b);
            else
    w(i,j)=w(i,j)+(16*q0/(pi^6*Dmn*m*n))*(sin(m*pi*x0/a)*sin(n*pi*y0/b)*sin(m*pi*c/(2*a))*sin(n*pi*d/(2*b)))*sin(m*pi*x/a)*sin(n*pi*y/b);
            end
            end 
        end
    end
end
x=linspace(0,a,nbpx);
y=linspace(0,b,nbpy);
w=w';
figure(1);
% Image de la déformée
surfc(x,y,w);
%shading interp
colorbar;
% figure(2);
% tri=delaunay(x,y);
% trimesh(tri,x,y,w);

% Calcul des contraintes maximale normales et hors plan

     

    case 5
        disp ('Plaque orthotrope appuyee en flambement bi-axiale')

% =========================================================================
% FLambement d'une plaque a  x   b   appuyée sur les 4 côtés
% =========================================================================

seq=[90,45,0,-45,-45,0,45,90];

nbp=length(seq);
ht=ep*length(seq);

k=0;

% initialisation des tableaux
A=zeros(3);
B=zeros(3);
D=zeros(3);
epsx=zeros(3,1);
sigl=zeros(3,nbp);
crit=zeros(1,nbp);
Nmax=zeros(1,nbp);

% Matrice de rigidité du pli dans le repère local
Slt=[1/El -nutl/Et 0 ; -nult/El 1/Et 0 ; 0 0 1/Glt];
Qlt=inv(Slt);

for i=1:(nbp+1)
    h(i)=-ht/2+(i-1)*ep;  
end

for i=1:nbp
    c=cos(seq(i)*pi/180);
    s=sin(seq(i)*pi/180);
    TS=[c*c s*s -2*s*c ; s*s c*c 2*s*c ; s*c -s*c c*c-s*s];
    TE=[c*c s*s -s*c ; s*s c*c s*c ; 2*s*c -2*s*c c*c-s*s]; % ou TE=(inv(TS))'     
% Matrice de rigidité/souplesse du pli dans le repère global

    Qxy=TS*Qlt*inv(TE);
    
    A=A+Qxy*(h(i+1)-h(i));
    B=B+(1/2)*Qxy*(h(i+1)^2-h(i)^2);
    D=D+(1/3)*Qxy*(h(i+1)^3-h(i)^3);
end

ABD=[A B; B D];
ABD

% Choix mini n=1
i=0;n=1;
for m=1:1;
    for a=25:5:800;
        for b=100:5:100;
    Ncr=(pi*pi/(a*a*(m*m+k*n*n*(a/b)^2)))*(D(1,1)*m^4+2*(D(1,2)+2*D(3,3))*m^2*n^2*(a/b)^2+D(2,2)*n^4*(a/b)^4);
    i=i+1;N11(1,i)=Ncr;
    xx1(1,i)=(D(2,2)/D(1,1))^(0.25)*a/b;
        end
    end
end
i=0;n=1;
for m=2:2;
    for a=25:5:800;
        for b=100:5:100;
    Ncr=pi*pi*(D(1,1)*m^4+2*(D(1,2)+2*D(3,3))*m^2*n^2*(a/b)^2+D(2,2)*n^4*(a/b)^4)/(a*a*(m*m+k*n*n*(a/b)^2));
    i=i+1;N21(1,i)=Ncr;
    xx2(1,i)=(D(2,2)/D(1,1))^(0.25)*a/b;
        end
    end
end
i=0;n=1;
for m=3:3;
    for a=25:5:800;
        for b=100:5:100;
    Ncr=pi*pi*(D(1,1)*m^4+2*(D(1,2)+2*D(3,3))*m^2*n^2*(a/b)^2+D(2,2)*n^4*(a/b)^4)/(a*a*(m*m+k*n*n*(a/b)^2));
    i=i+1;N31(1,i)=Ncr;
    xx3(1,i)=(D(2,2)/D(1,1))^(0.25)*a/b;
        end
    end
end
i=0;n=1;
for m=4:4;
    for a=25:5:800;
        for b=100:5:100;
    Ncr=pi*pi*(D(1,1)*m^4+2*(D(1,2)+2*D(3,3))*m^2*n^2*(a/b)^2+D(2,2)*n^4*(a/b)^4)/(a*a*(m*m+k*n*n*(a/b)^2));
    i=i+1;N41(1,i)=Ncr;
    xx4(1,i)=(D(2,2)/D(1,1))^(0.25)*a/b;
        end
    end
end
for m=5:5;
    for a=25:5:800;
        for b=100:5:100;
    Ncr=pi*pi*(D(1,1)*m^4+2*(D(1,2)+2*D(3,3))*m^2*n^2*(a/b)^2+D(2,2)*n^4*(a/b)^4)/(a*a*(m*m+k*n*n*(a/b)^2));
    i=i+1;N51(1,i)=Ncr;
    xx5(1,i)=(D(2,2)/D(1,1))^(0.25)*a/b;
        end
    end
end
figure;
plot(xx1,N11,xx2,N21,xx3,N31,xx4,N41,xx5,N51);
xlim([0 6])
xlabel('(D22/D11)^0^.^2^5.R');
ylim([100 300])
ylabel('Effort Critique (N/mm)');

grid on;


    case 6
        
        disp ('Poutre stratifiée ou Sandwich en flexion');
        
        disp('pas encore fait!');
        

end






