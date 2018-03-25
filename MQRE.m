% clear all; clc; close all;

%% Modelo AR com minimos quadrados recurssivos
 nit=length(y); % numero de amostras
 t=0:Ts:nit*Ts-Ts; % Time vector based on Ts.
%% 
% close all;
N=6; % tamanho 
Y=y';
U=u';
e(1:5)=0;
theta_hat= zeros(N,1);  % valor inicial do vetor de parâmetros. 
P=1000*eye(N,N);       % valor inicial da matriz de covariância. 
% nit=size(u,1);               % numero de amostras 
                                    % entrada-> vetor "u" do SBPA 
                                    % saída-> vetor "y" obtido com SBPA 
Yhat=zeros(size(y));    % vetor de saídas estimadas 
PHI=zeros(N,1);           % vetor de dados 
THETA=zeros(size(u,1),N); %vetor de parâmetros estimados 
ek=0;                            % erro de previsão 
W(1:5)=0;
yhat=0;                         % saida estimada atual 
K=zeros(6,1);               % ganho do estimador 
TRACOP=zeros(size(u,1),1);
lambda=1; % fator de esquecimento
I=eye(N); % Inicializar a matriz identidade
for t=5:nit;  
 % (1) :  "mede" os valores de entrada e saída da planta;    
         % feito através do SBPA 
 % (2) : Forma o vetor de dados PHI;   
     
  PHI = [-Y(t-1) -Y(t-2) U(t-1) U(t-2) W(t-1) W(t-2)]';  
 % (3) : Calcula a saída estimada pelo modelo no instante k, yhat;    
   
  yhat= PHI' * theta_hat; 
   
  %guarda as saídas obtidas no vetor Yhat 
  Yhat(t)=yhat; 
   
  % (4) : Calcula o valor do erro de estimação no instante k; 
  
  ek=y(t)-yhat;   e(t)=ek;
  
  % (5) :  Calcula o vetor de ganho do estimador:  
  
  denominador= lambda + PHI' * P * PHI; 
  K=(1/denominador)*P*PHI; 
   
  % (6): atualiza estimativa do vetor de parâmetros theta_hat; 
  
  theta_hat= theta_hat + K*ek;  
   
 % ruído brando, processo estocástico 
  W(t)=y(t)-(PHI' * theta_hat);
  
  THETA(t,1)=theta_hat(1);  %armazena trajetória de a1_h%at.  
  THETA(t,2)=theta_hat(2);  %armazena trajetória de a2_hat.  
  THETA(t,3)=theta_hat(3);  %armazena trajetória de a3_hat.  
  THETA(t,4)=theta_hat(4);  %armazena trajetória de b0_hat.  
  THETA(t,5)=theta_hat(5);  %armazena trajetória de b1_hat.  
  THETA(t,6)=theta_hat(6);  %armazena trajetória de b2_hat.  
  
  % passo 7: atualiza a matriz de covariância; 
   
  P = (I-K*PHI')*P/lambda;
  TRACOP(t)= trace(P); 
   
end; %final do laço for principal.  
% plot(THETA)     %parâmetros estimados
a1=theta_hat(1);a2=theta_hat(2);b0=theta_hat(3);b1=theta_hat(4);
c1=theta_hat(5); c2=theta_hat(6); 

clear ys
ys(1:5)=0;
es(1:10)=0;
for k=5:nit
    
    ys(k)=-a1*y(k-1)-a2*y(k-2)+ b0*u(k-1) + b1*u(k-2) +e(k)+c1*e(k-1)+c2*e(k-2) ;  
    es(k)=y(k)-ys(k);
end
% figure;
plot(y,'b'); hold
plot(ys,'r');
var(W)  
%% Erro Quadrático 
e_AR = y'-ys'; % Erro da saída real e a saída estimada
e2_AR = e_AR.^2; % Erro quadrático
SEQ_AR = sum(e2_AR); % Somatório dos erros quadráticos

%% Coeficiente de Correlação Múltipla RMS
m1_AR = sum((y-mean(ys)).^2); % Soma a diferença entre cada amostra e a média das amostras elevadas ao quadrado armazenadas
R2_AR = (1-(SEQ_AR/m1_AR)) *100 % Coeficiente de correlação múltipla
title(['fit de correlação= ' num2str(R2_AR)]); hold;
%% Forma canonica observavel em espaço de estados 
% Parameters in State Space
A = [1 a1 a2]
B = [b0 b1 0]
C = [1 c1 c2 ]
Gz=tf(B,A,Ts)