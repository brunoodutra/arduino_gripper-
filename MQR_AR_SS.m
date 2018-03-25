% clear all; close all; clc;
function [A,C]=MQR_AR_SS(y,Ts);
% Importa os Dados
% clear all; clc; close all;
% [X] = importdata('teste3.TXT');
% [t, y1, y2, y3, y4] = textread('teste1.TXT', ...
% '%f %f %f %f %f',60001);
% Ts = mean(diff(t));
% Win=80;
% y1=abs(y1);
% y1=MAV(Win,length(y1),y1);
% y=y1;
u=y*0;% modelo autoregressivo
%% Modelo AR com minimos quadrados recurssivos
 nit=length(y); % numero de amostras
 t=0:Ts:nit*Ts-Ts; % Time vector based on Ts.
%% 
% close all;
N=8; % tamanho 
Y=y';
U=y'*0;
e(1:5)=0;
theta_hat= zeros(N,1);  % valor inicial do vetor de par�metros. 
P=1000*eye(N,N);       % valor inicial da matriz de covari�ncia. 
% nit=size(u,1);               % numero de amostras 
                                    % entrada-> vetor "u" do SBPA 
                                    % sa�da-> vetor "y" obtido com SBPA 
Yhat=zeros(size(y));    % vetor de sa�das estimadas 
PHI=zeros(N,1);           % vetor de dados 
THETA=zeros(size(u,1),N); %vetor de par�metros estimados 
ek=0;                            % erro de previs�o 
W(1:5)=0;
yhat=0;                         % saida estimada atual 
K=zeros(6,1);               % ganho do estimador 
TRACOP=zeros(size(u,1),1);
lambda=1; % fator de esquecimento
I=eye(N); % Inicializar a matriz identidade
for t=5:nit;  
 % (1) :  "mede" os valores de entrada e sa�da da planta;    
         % feito atrav�s do SBPA 
 % (2) : Forma o vetor de dados PHI;   
     
  PHI = [-Y(t-1) -Y(t-2) -Y(t-3) -Y(t-4) U(t-1) U(t-2) U(t-3) U(t-4)]';  
 % (3) : Calcula a sa�da estimada pelo modelo no instante k, yhat;    
   
  yhat= PHI' * theta_hat; 
   
  %guarda as sa�das obtidas no vetor Yhat 
  Yhat(t)=yhat; 
   
  % (4) : Calcula o valor do erro de estima��o no instante k; 
  
  ek=y(t)-yhat;   e(t)=ek;
  
  % (5) :  Calcula o vetor de ganho do estimador:  
  
  denominador= lambda + PHI' * P * PHI; 
  K=(1/denominador)*P*PHI; 
   
  % (6): atualiza estimativa do vetor de par�metros theta_hat; 
  
  theta_hat= theta_hat + K*ek;  
   
 % ru�do brando, processo estoc�stico 
  W(t)=y(t)-(PHI' * theta_hat);
  
  THETA(t,1)=theta_hat(1);  %armazena trajet�ria de a1_h%at.  
  THETA(t,2)=theta_hat(2);  %armazena trajet�ria de a2_hat.  
  THETA(t,3)=theta_hat(3);  %armazena trajet�ria de a3_hat.  
  THETA(t,4)=theta_hat(4);  %armazena trajet�ria de b0_hat.  
  THETA(t,5)=theta_hat(5);  %armazena trajet�ria de b1_hat.  
  THETA(t,6)=theta_hat(6);  %armazena trajet�ria de b2_hat.  
  THETA(t,7)=theta_hat(7);  %armazena trajet�ria de c1_hat.  
  THETA(t,8)=theta_hat(8);  %armazena trajet�ria de c2_hat.  
%   THETA(t,9)=theta_hat(9);  %armazena trajet�ria de c3_hat.  
%    
  % passo 7: atualiza a matriz de covari�ncia; 
   
  P = (I-K*PHI')*P/lambda;
  TRACOP(t)= trace(P); 
   
end; %final do la�o for principal.  
% plot(THETA)     %par�metros estimados
a1=theta_hat(1);a2=theta_hat(2);a3=theta_hat(3);a4=theta_hat(4);
c1=theta_hat(5); c2=theta_hat(6); c3=theta_hat(7); c4=theta_hat(8)
% c1=theta_hat(7); c2=theta_hat(8); c3=theta_hat(9);
clear ys
ys(1:5)=0;
es(1:10)=0;
for k=5:nit
    
    ys(k)=-a1*y(k-1)-a2*y(k-2)-a3*y(k-3)-a4*y(k-4)+e(k)+c1*e(k-1)+c2*e(k-2)+c3*e(k-3)+c4*e(k-4) ;  
    es(k)=y(k)-ys(k);
end
% figure;
plot(y,'b'); hold
plot(ys,'r');
var(W)  
%% Erro Quadr�tico 
e_AR = y'-ys'; % Erro da sa�da real e a sa�da estimada
e2_AR = e_AR.^2; % Erro quadr�tico
SEQ_AR = sum(e2_AR); % Somat�rio dos erros quadr�ticos

%% Coeficiente de Correla��o M�ltipla RMS
m1_AR = sum((y-mean(ys)).^2); % Soma a diferen�a entre cada amostra e a m�dia das amostras elevadas ao quadrado armazenadas
R2_AR = (1-(SEQ_AR/m1_AR)) *100 % Coeficiente de correla��o m�ltipla
title(['fit de correla��o= ' num2str(R2_AR)]); hold;
%% Forma canonica observavel em espa�o de estados 
A=[1 a1 a2 a3 a4];
dA = conv([1 -1],A);
na = length(A)-1;

% Parameters in State Space
A = [-A(2) 1 0 0;
     -A(3) 0 1 0;
     -A(4) 0 0 1
     -A(5) 0 0 0];

% B = [b0 b1 b2 0]';
G = -dA(2:end)';
C = [1 0 0 0];
x(:,:,1)=[0,0,0,0]';
Qfk = diag([1 1 1 1]); Rfk =1;

  G = dlqr(A',[1 c1 c2 c3]',Qfk,Rfk)';
% G = [(c1-dA(2)) (c2-dA(3)) (c2-dA(4)) (c3-dA(5))]';% Gama 
  for k=2:nit,
      x(:,:,k) =(A-G*C)*x(:,:,k-1) + G*y(k-1);
%       es(k)=y(k)-ys(k);
  end
%   plot(x(1,:),'g')

% theta_hat'
end