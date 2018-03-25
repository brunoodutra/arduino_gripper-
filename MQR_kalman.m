% % clear all; close all; clc;
function [A,C]=MQR_kalman(y,Ts);
% Importa os Dados
% [X] = importdata('teste3.TXT');
% [t, y1, y2, y3, y4] = textread('teste1.TXT', ...
% '%f %f %f %f %f',60001);
% 
%  y=y1;

u=y*0;
%% Modelo AR com minimos quadrados
 nit=length(y); % numero de amostras
 t=0:Ts:nit*Ts-Ts; % Time vector based on Ts.
%% 
% close all;
N=6; % tamanho 
Y=y';
U=y'*0;
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
for t=4:nit;  
 % (1) :  "mede" os valores de entrada e saída da planta;    
         % feito através do SBPA 
 % (2) : Forma o vetor de dados PHI;   
     
  PHI = [-Y(t-1) -Y(t-2) -Y(t-3)  W(t-1) W(t-2) W(t-3)]';  
 % (3) : Calcula a saída estimada pelo modelo no instante k, yhat;   
   
  yhat= PHI' * theta_hat; 
   
  %guarda as saídas obtidas no vetor Yhat 
  Yhat(t)=yhat; 
   
  % (4) : Calcula o valor do erro de estimação no instante k; 
   
  ek=y(t)-yhat;
  e(t)=ek;
  % (5) :  Calcula o vetor de ganho do estimador: 
   
  denominador= lambda + PHI' * P * PHI; 
  K=(1/denominador)*P*PHI; 
   
  % passo 6: atualiza estimativa do vetor de parâmetros theta_hat; 
   
  theta_hat= theta_hat + K*ek; 
   
  W(t)=y(t)-(PHI' * theta_hat);
  
  THETA(t,1)=theta_hat(1);  %armazena trajetória de a1_hat.  
  THETA(t,2)=theta_hat(2);  %armazena trajetória de a2_hat.  
  THETA(t,3)=theta_hat(3);  %armazena trajetória de a3_hat.  
  THETA(t,4)=theta_hat(4);  %armazena trajetória de b0_hat.  
  THETA(t,5)=theta_hat(5);  %armazena trajetória de b1_hat.  
  THETA(t,6)=theta_hat(6);  %armazena trajetória de b2_hat.  
%   THETA(t,7)=theta_hat(7);  %armazena trajetória de c1_hat.  
%   THETA(t,8)=theta_hat(8);  %armazena trajetória de c2_hat.  
%   THETA(t,9)=theta_hat(9);  %armazena trajetória de c3_hat.  
%    
  % passo 7: atualiza a matriz de covariância; 
   
%   P= P - K*(P*PHI)'; 
  P = (I-K*PHI')*P/lambda;
  TRACOP(t)= trace(P); 
   
end; %final do laço for principal.  
% plot(THETA)     %parâmetros estimados
a1=theta_hat(1);a2=theta_hat(2);a3=theta_hat(3);
c1=theta_hat(4); c2=theta_hat(5); c3=theta_hat(6);
% c1=theta_hat(7); c2=theta_hat(8); c3=theta_hat(9);
clear ys
ys(1:4)=0;
es(1:10)=0;
for k=4:nit
    
    ys(k)=-a1*ys(k-1)-a2*ys(k-2)-a3*ys(k-3)+e(k)+c1*e(k-1)+c2*e(k-2)+c3*e(k-3);    
    es(k)=y(k)-ys(k);
end
% figure;
% plot(y,'b'); hold
% plot(ys,'r');

%% Forma canonica observavel em espaço de estados 
A=[1 a1  a2 a3];
dA = conv([1 -1],A);
% B=[0 b0 b1 b2];
C=[1  c1  c2 c3];
na = length(A)-1;

% Parameters in State Space
A = [-dA(2) 1 0 0;
     -dA(3) 0 1 0;
     -dA(4) 0 0 1
     -dA(5) 0 0 0];

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
  
%    plot(x(1,:))

theta_hat'
end