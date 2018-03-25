clear all;clc; close all;

s=tf('s');


[x t]=step(1/(2*s+1),15);


% x=sin(t)
Ts = mean(diff(t));
nit=length(x)
%% ruído 1
variancia = 0.05;
z = rand(nit,1); z = 2*z-1; % Variável auxiliar
xi = sqrt(variancia*3)*z; % Ruído Branco
%% ruído2
z1 = rand(nit,1); z1 = 2*z1-1; % Variável auxiliar
xi1 = sqrt(variancia*3)*z1; % Ruído Branco

y1=x+xi;
y2=x+xi1;

[Ad,Cd]=MQR_AR_SS(x,Ts);
% [Ad_1,Cd_1]=MQR_AR_SS(y2,Ts); 
% Ad=[Ad(:,1) Ad_1(:,1) eye(4,2)]

Temp=t;
nit = length(y1);


%% Filtro Kalman dual ao LQR de horizonte infinito 
Qfk1 = diag([100.0000 1.0000 1.0000 1.0000]); Rfk1 =5e+4;
%% Calculo recursivo da matriz de covariancia P e do ganho L do filtro de kalman     
p1(:,:,1)=100*diag([1 1 1 1]);
for k=2:300;
    p1(:,:,k)= Ad*p1(:,:,k-1)*Ad' -Ad*p1(:,:,k-1)*Cd'*inv(Cd*p1(:,:,k-1)*Cd'+Rfk1)*Cd*p1(:,:,k-1)*Ad'+Qfk1;
    P1=p1(:,:,k);
end  
L1=(Ad*P1*Cd'*inv(Cd*P1*Cd'+Rfk1));
% % % % % % % % % % % % % % % % % % % % % % % % 
%% fusão dos dados y1 e y2 referêntes ao musculo flexor
tempsimul=nit;  
R=1;
% C=(Ad-L*Cd)
I=diag([1 1 1 1]);
xa=[0;0;0;0];Pa_=0; Pa=P1; yhat(1)=0;
xi=[0;0;0;0]; Pi_=0; Pi=P1; y(1)=0; 
P1_=Pa+Pi; P11=P1_; x1=[0;0;0;0]; M1=1000000; M2=1000000; m1=xa; m2=xi;
emg1(1)=0;emg2(1)=0;
for i=2:tempsimul
   emg1(i)=y1(i);
   
   if i>30 && i<34
%         emg2(i)=y2(i)*0;
        emg2(i)=y2(i);

    else
        emg2(i)=y2(i);
   end

   
% Filtro de kalman 1
xi=(Ad-L1*Cd)*xi + L1*emg1(i);
y(i)= Cd*xi;

% Filtro de kalman 2
xa=(Ad-L1*Cd)*xa + L1*emg2(i);
yhat(i)= Cd*xa;

%% Fuzão  
% Covariança do fusor
%     P1_=P1_+(Pa-M1)+(Pi-M2)
    x1= inv(P1_)*(P11*x1+(Pa*xa - M1*m1)+(Pi*xi - M2*m2)) ;
    f1(i)=Cd*x1;
    M1=Pa; M2=Pi; m1=xa; m2=xi; P11=P1_
%            
% x1= inv(P1_)*( P11*x1 + Pa*(xa-m1)+Pi*(xi-m2) );
%     
%     f1(i)=Cd*x1;
%     M1=Pa; M2=Pi; m1=xa; m2=xi; P11=P1_
end

%% [R,P]=corrcoef([emg1' emg2' y1'])
figure
subplot(211)
plot(t,emg1,'k',t,emg2,'g'); hold;
plot(t,x,'b','linewidth',2);
legend('Sinal medido 1','Sinal medido 2','Sinal Original' );
ylabel('Volts');
xlabel('Tempo (s)');
title('Teste em ambiente controlado')
subplot(212)
plot(t,x,'b','linewidth',2);hold
plot(t,y,'k',t,yhat,'g','linewidth',1.5); 
plot(t,f1,'r','linewidth',2);
legend('Sinal Original','Sinal_F_k1','Sinal_F_k2','Sinal_F_D_K' );
ylabel('Volts');
xlabel('Tempo (s)');

[R1,P1]=corrcoef([x y']);
[R2,P2]=corrcoef([x yhat']);
[R3,P3]=corrcoef([x f1']);
indices=[R1(2) R2(2) R3(2)]

%% Erro Quadrático yest1
e_1 = x-y'; % Erro da saída real e a saída estimada
e2_1 = e_1.^2; % Erro quadrático
SEQ_1 = sum(e2_1); % Somatório dos erros quadráticos
%% Coeficiente de Correlação Múltipla RMS
m1_1 = sum((x-mean(y)).^2); % Soma a diferença entre cada amostra e a média das amostras elevadas ao quadrado armazenadas
R2_1 = (1-(SEQ_1/m1_1)) *100; % Coeficiente de correlação múltipla
R2_1;

%% Erro Quadrático yest2
e_2 = x-yhat'; % Erro da saída real e a saída estimada
e2_2 = e_2.^2; % Erro quadrático
SEQ_2 = sum(e2_2); % Somatório dos erros quadráticos
%% Coeficiente de Correlação Múltipla RMS
m1_2 = sum((x-mean(yhat)).^2); % Soma a diferença entre cada amostra e a média das amostras elevadas ao quadrado armazenadas
R2_2 = (1-(SEQ_2/m1_2)) *100; % Coeficiente de correlação múltipla
R2_2;

%% Erro Quadrático yest3
e_3 = x-f1'; % Erro da saída real e a saída estimada
e2_3 = e_3.^2; % Erro quadrático
SEQ_3 = sum(e2_3); % Somatório dos erros quadráticos
%% Coeficiente de Correlação Múltipla RMS
m1_3 = sum((x-mean(f1)).^2); % Soma a diferença entre cada amostra e a média das amostras elevadas ao quadrado armazenadas
R2_3 = (1-(SEQ_3/m1_3)) *100; % Coeficiente de correlação múltipla
R2_3;

indices2=[R2_1 R2_2 R2_3]