%#########################################################################%
%---------------TRAITEMENT D'UN SIGNAL ELECTROGARDIOGRAMME----------------%
%---------------------------EL ASRI YASSINE-------------------------------%
%-----------------------MASTER GE A TEMPS AMENAGE-------------------------%
%---------------------MODULE DE TRAITEMENT DE SIGNAL----------------------%
%---------------------SUPERVISER PAR Mr.ELMHAMDI--------------------------%
%#########################################################################%

clc
clear all
close all
load ('I03.mat')
s=val(1,:);                                 %Choix du signal N°1 
whitebg('k');                               %Affichage des figures en noir
figure(1), plot(s),title('signal d entrée') %Affichage du signal d'entrée
M=length(s)                                 %Affichage de la languer du signal
%--------------------------------------------------------------------------

%------------------------TRASFORMEE DE FOURIER RAPIDE----------------------
F=fft(s)                               %F c'est La FFT du signal d'entrée
figure(2), plot(abs(F)); title('FFT'); %Affichage du FFT signal d'entrée
%--------------------------------------------------------------------------

%---------------------------Création du bruit------------------------------
bruit=wgn(1,M,1);                     %Creation d'un bruit blanc gaussien          
sb=bruit;                             %Le signal de Bruit
figure(3), plot(sb);title('BRUIT');   %Affichage du FFT signal d'entrée
%--------------------------------------------------------------------------

%--------------------------Filtrage du signal------------------------------
fe=1000;                                 %Frequence d'échantillonnage
Te=1/fe                                  %Periode d'echantillonnage
fc=50;                                   %Frequence de coupure
O =20;                                   %Ordre du filtre
fcn=fc/(fe/2);                           %Frequence de coupure
[num,den]=butter(15,fcn);                %Fitrage
Hb=dsp.IIRFilter('Numerator',num,'Denominator',den);
d = step(Hb,s');                         %d est le signal filtré desiré
figure(4);plot(d);title('Signal desiré');%Affichage du signal filtré
figure(5);plot(d,'ro'); hold on; plot (d);title('Signal desiré');    
%--------------------------------------------------------------------------

%---------------------------Filtre Adaptatif LMS--------------------------- 
delta=0.00005; % Facteur d'adaptation
hlms=dsp.LMSFilter(O,'StepSize',delta); %POUR REDUIRE L'ERREUR IL FAUT JOUER SUR DELTA
[y,e,w] = step(hlms,sb',d);
figure(6),plot(1:M,[d,y,e]),title('Transposition des Signaux desiré');
legend('Signal desiré','Sortie','Erreur');
%--------------------------------------------------------------------------

%---------------------Erreur prédite et la sortie prédite------------------
a=lpc(s,O);           %Linear prediction filter coefficients a l'ordre O
err=filter(a,1,sb);   %extraction du Signal predit
rest=filter(1,a,err); %extraction du Bruit predit
figure(7),plot(1:M-1,rest(1:M-1)),title('Signal de bruit predit');%Affichage du signal predit
figure(8),plot(1:M-1,err(1:M-1)),title('Erreur predit'); %Affichage de l'erreur predit
%--------------------------------------------------------------------------

%-------------------------Evolution des coefficients-----------------------
N=50;                           
h=zeros(N+1,O);
e(O)=0;
for i=O+1:N;
 y(i)=0;
 for j=1:20;
 y(i)=y(i)+h(i,j)*s(1,i-j+1);
 h(i+1,j)=h(i,j)+delta*e(i-1)*s(1,i-j+1);
 end;
end;
figure(8);
plot(0:N-O,h(20:N,1),'r'), hold on;
plot(0:N-O,h(20:N,2),'r'), hold on;
plot(0:N-O,h(20:N,3),'r'), hold on;
plot(0:N-O,h(20:N,4),'r'), hold on;
plot(0:N-O,h(20:N,5),'r'), hold on;
plot(0:N-O,h(20:N,6),'r'), hold on;
plot(0:N-O,h(20:N,7),'g'), hold on;
plot(0:N-O,h(20:N,8),'g'), hold on;
plot(0:N-O,h(20:N,9),'g'), hold on;
plot(0:N-O,h(20:N,10),'g'), hold on;
plot(0:N-O,h(20:N,11),'g'), hold on;
plot(0:N-O,h(20:N,12),'b'), hold on;
plot(0:N-O,h(20:N,13),'b'), hold on;
plot(0:N-O,h(20:N,14),'b'), hold on;
plot(0:N-O,h(20:N,15),'b'), hold on;
plot(0:N-O,h(20:N,16),'b'), hold on;
plot(0:N-O,h(20:N,17),'b'), hold on;
plot(0:N-O,h(20:N,18),'b'), hold on;
plot(0:N-O,h(20:N,19),'b'), hold on;
plot(0:N-O,h(20:N,20),'b'), hold off;
grid;
title('Evolution des coefficients ');
%--------------------------------------------------------------------------

%-----------------------DETECTION DES R PEAKS-----------------------
Compteur_de_battements = 0;
for k = 2 : M-1
    if(s(k) > s(k-1) & s(k) > s(k+1) & s(k)> 4000)
        k
        disp('Un Peak R Domminant Trouvé');
        Compteur_de_battements = Compteur_de_battements + 1
        fs = (s(k)-s(k-1))
    end
end
%--------------------------------------------------------------------------

%-------------------Calcul de nombre de battement par minute---------------
fs=mean(fs);                %Calcul de la valeur moyenne de la frequence
durer_en_secondes = M/fs;   %Calcul de la durée du signal en secondes 
durer_en_minutes = durer_en_secondes/60;%Calcul de la durée du signal en Min 
Moyenne_des_Battements = Compteur_de_battements/durer_en_minutes %Calcul de nombre de battement en secondes 
%--------------------------------------------------------------------------

%-----------------------Affichage des echantillons-------------------------
S_COEUR = d(2320:2465);
figure(9), plot(S_COEUR);title('SIGNAL_COEUR');
P = d(2320:2362);
Fp = 1/(2362-2320)
figure(10), plot(P); title('SIGNAL_P');
QRS = d(2370:2412);
Fqrs = 1/(2412-2370)
figure(11), plot(QRS); title('SIGNAL_QRS');
T = d(2440:2466);
Ft = 1/(2466-2440)
figure(12), plot(T); title('SIGNAL_T');
%--------------------------------------------------------------------------

