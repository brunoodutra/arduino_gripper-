clear all;clc; close all;
format long
% [X] = importdata('teste3.TXT');
load('coleta6.mat')
% [t, y1, y2, y3, y4] = textread('teste1.TXT', ...
% '%f %f %f %f %f',60001);
% analise_dos_dados(t,y4);
%%
% % Ts = mean(diff(t)); 
% % Ts=5.0000e-04; % tempo de amostragem
% % Ts=5.0000
% % Y1=y1;
% % y2(1:500)=y2(501:1000);
% nit=length(y2);
% Win=1;
% 
% y1=abs(y1);
% y1=MAV(Win,nit,y1);
y1=Norm(y1,5,0);
% % y1 = detrend(y1);
% 
% y2=abs(y2);
% y2=MAV(Win,nit,y2);
y2=Norm(y2,5,0);
% % y2 = detrend(y2);
% 
% y3=abs(y3);
% y3=MAV(Win,nit,y3);
y3=Norm(y3,5,0);
% % y3 = detrend(y3);
% 
% y4=abs(y4);
% y4=MAV(Win,nit,y4);
y4=Norm(y4,5,0);
% % y4 = detrend(y4);
% 
% Ts=Win*Ts
%% modelo discretizado 
% Ad    =      [  0.98717,     -0.11654,     0.047194,   -0.0039135;
%               -0.077789,     0.042169,      0.73997,      0.43509;   
%               -0.027153,     -0.57941,     -0.26153,      0.73948;
%               -0.0027486,      0.30031,     -0.61597,      0.22622];
% 
% Cd =[   13757,       1466.3,       1148.2,       1444.2];


[Ad,Cd]=MQR_AR_SS(y1,Ts);
[Ad_1,Cd_1]=MQR_AR_SS(y2,Ts); 
% Ad=[Ad(:,1) Ad_1(:,1) eye(4,2)]


[Ad2,Cd2]=MQR_AR_SS(y3,Ts);
[Ad2_1,Cd2_1]=MQR_AR_SS(y4,Ts);
% Ad2=[Ad2(:,1) Ad2_1(:,1) eye(4,2)]

t =0:Ts:length(y1)*Ts-Ts;
Temp=t;
nit = length(y1);


%% Filtro Kalman dual ao LQR de horizonte infinito 
Qfk1 = diag([020.0000 1.0000 1.0000 1.0000]); Rfk1 =5e+2;
Qfk2 = diag([020.0000 1.0000 1.0000 1.0000]); Rfk2 =5e+2;
% Qfk1 = diag([010.0000 1.0000 1.0000 1.0000]); Rfk1 =1e+3;
% Qfk2 = diag([010.0000 1.0000 1.0000 1.0000]); Rfk2 =1e+3;
%% Calculo recursivo da matriz de covariancia P e do ganho L do filtro de kalman     
p1(:,:,1)=100*diag([1 1 1 1]);
p2(:,:,1)=100*diag([1 1 1 1]);
for k=2:300;
    p1(:,:,k)= Ad*p1(:,:,k-1)*Ad' -Ad*p1(:,:,k-1)*Cd'*inv(Cd*p1(:,:,k-1)*Cd'+Rfk1)*Cd*p1(:,:,k-1)*Ad'+Qfk1;
    P1=p1(:,:,k);
    
    p2(:,:,k)= Ad2*p2(:,:,k-1)*Ad2' -Ad2*p2(:,:,k-1)*Cd'*inv(Cd*p2(:,:,k-1)*Cd'+Rfk2)*Cd*p2(:,:,k-1)*Ad2'+Qfk2;
    P2=p2(:,:,k);
end  
L1=(Ad*P1*Cd'*inv(Cd*P1*Cd'+Rfk1));
L2=(Ad2*P2*Cd'*inv(Cd*P2*Cd'+Rfk2));
% % % % % % % % % % % % % % % % % % % % % % % % 

%% análise de robustez filtro de kalman 

%% Filtro Kalman

% Analise de robustez do FK
    [num,den] = ss2tf(Ad-L1*Cd, L1, Cd,0);
    Tfk = tf(num,den,Ts);
    
    ganhos=dsigma(Tfk,Ts);
    g=ganhos;
    figure
    dsigma(Tfk,Ts)
    Mt = max(ganhos); clear ganhos;

    Sfk = tf(1,den,Ts);
          ganhos=dsigma(Sfk,Ts);
          Ms = max(ganhos); clear ganhos;
         
         disp('FK:');
         % MG
           MGfk = max( (Ms/(Ms-1)), (1+ (1/Mt)) );
           MGfk_dB = 20*log10(MGfk)
         % MF
           MFfk = max( (2*asin(1/(2*Ms))*(180/pi)), ...
               (2*asin(1/(2*Mt))*(180/pi)) )
       
     figure
      dsigma(Sfk,Ts,'inv');
      hold on;
      dsigma(Tfk,Ts);
      title(['{Filtro de Kalman: MG >= }' num2str(MGfk_dB) 'dB' ...
          '{ e MF >= }' num2str(MFfk) '{\circ}.']);
          legend('| T |');
         ylabel('Ganhos Principais');
     
           
           
     
%     % Sensibilidade complementar: GP_T : ganhos principais de T
%       [GP_T,w] = dsigma(Sfk,Ts); GP_T=mag2db(GP_T);
%       Mt = max( max(GP_T) ); % "Amplitude ratio" de T
%       MG_T = 1+(1/Mt);
%       MF_T = 2*asin( 1/(2*Mt) )*(180/pi);
%       
%     % Sensibilidade: GP_S : ganhos principais de S
%       [GP_S,w] = dsigma(Sfk,Ts,'inv');GP_S=mag2db(GP_S);
%       Ms = max( max(GP_S(1))); % "Amplitude ratio" de S
%       MG_S = Ms/(Ms-1);
%       MF_S = 2*asin( 1/(2*Ms) )*(180/pi);
%       
%     % Resultados MG e MF pelas curvas de sensibilidade
%       MG_fk = max(MG_T,MG_S);
%       MF_fk = max(MF_T,MF_S);
%       figure(10);
%       dsigma(Sfk,Ts);hold;dsigma(Sfk,Ts,'inv');
%       title(['{Filtro de Kalman: MG >= }' num2str(MG_fk) 'dB' ...
%           '{ e MF >= }' num2str(MF_fk) '{\circ}.']);
%       legend('|T|','|S|');
%       ylabel('Ganhos Principais');
% 

%%
%% ruído
variancia = 0.005;
x = rand(nit,1); x = 2*x-1; % Variável auxiliar
ri = sqrt(variancia*5)*x; % Ruído Branco
ri=ri;
%% fusão dos dados y1 e y2 referêntes ao musculo flexor
tempsimul=nit;  
R=1;
% C=(Ad-L*Cd)
I=diag([1 1 1 1]);
xa=[0;0;0;0];Pa_=0; Pa=P1; y_est1(1)=0;
xi=[0;0;0;0]; Pi_=0; Pi=P1; y_est2(1)=0;
P1_=Pa+Pi; P11=P1;x1=[0;0;0;0]; M1=1000; M2=1000; m1=xa; m2=xi;
emg1(1)=0;emg2(1)=0;
for i=2:tempsimul
   emg1(i)=y1(i);
   
   if i>1740 && i<1780
        emg2(i)=0.1+ri(i);
%         emg2(i)=y2(i);

    else
        emg2(i)=y2(i);
   end

   
% Filtro de kalman 1
xi=(Ad-L1*Cd)*xi + L1*emg1(i);
y_est1(i)= Cd*xi;

% Filtro de kalman 2
xa=(Ad-L1*Cd)*xa + L1*emg2(i);
y_est2(i)= Cd*xa;

%% Fuzão  
% Covariança do fusor
%     P1_=P1_+(Pa-M1) + (Pi-M2)
    x1= inv(P1_)*( P1_*x1 + (Pa*xa - M1*m1)+(Pi*xi - M2*m2) );
    y1(i)=Cd*x1;
    M1=Pa; M2=Pi; m1=xa; m2=xi; P11=P1_;
           
% x1= inv(P1_)*( P1_*x1 + Pa*(xa-m1)+Pi*(xi-m2) );
%     
%     y1(i)=Cd*x1;
%     M1=Pa; M2=Pi; m1=xa; m2=xi;
end
%% fusão dos dados y3 e y4 referêntes ao musculo extensor
tempsimul=nit;  
% C=(Ad-L*Cd)
I=diag([1 1 1 1]);
xa=[0;0;0;0];Pa_=0; Pa=P2; y_est3(1)=0;
xi=[0;0;0;0]; Pi_=0; Pi=P2; y_est4(1)=0;
P2_=Pa+Pi; x1=[0;0;0;0]; M1=1000; M2=1000; m1=xa; m2=xi;
emg3(1)=0;emg4(1)=0;
for i=2:tempsimul
   emg3(i)=y3(i);
     if i>800 && i<=850
        emg4(i)=0.1 +ri(i);
    else
        emg4(i)=y4(i);
   end
   
    
xi=(Ad2-L2*Cd2)*xi + L2*emg3(i);
y_est3(i)= Cd2*xi;
     
xa=(Ad2-L2*Cd2)*xa + L2*emg4(i);
y_est4(i)= Cd2*xa;

%% fuzão
    % 2. Covariança do erro a frente
%     P1_=P1_+(Pa-M1) + (Pi-M2);
%     x1= inv(P1_)*( P1_*x1 + (Pa*xa - M1*m1)+(Pi*xi - M2*m2) );
%     y2(i)=Cd*x1;
%     M1=Pa; M2=Pi; m1=xa; m2=xi;

x1= inv(P2_)*( P2_*x1 + Pa*(xa-m1)+Pi*(xi-m2) );
    
    y2(i)=Cd*x1;
    M1=Pa; M2=Pi; m1=xa; m2=xi;

end

%%
    % figure
    % t1=emg1(680:1000);
    % t2=y_est1(680:1000);
    % t3=0:Ts:Ts*length(t2)-Ts;
    % plot(t3,t1,'b'); hold;
    % plot(t3,t2,'r','linewidth',2);
    % legend('Sinal_r_e_a_l','Sinal_F_K');
    % ylabel('Amplitude Normalizada');
    % xlabel('Tempo (s)');
    % title('Sinal EMG do flexor do punho vs. Sinal filtrado por FK');
%%
[R,P]=corrcoef([emg1' emg2' y1])
figure
subplot(221)
% plot(t,emg1,'b'); hold;
plot(t,emg1,'b'); hold;
legend('flexão canal1');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);
subplot(222)
% plot(t,emg2,'g'); hold;
plot(t,emg2,'b');hold
legend('flexão canal2');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);
subplot(223)
% plot(t,emg3,'b'); hold
plot(t,emg3,'b');hold
legend('Extensão canal3');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);
subplot(224)
% plot(t,emg4,'g'); hold
plot(t,emg4,'b');hold
legend('Extensão canal4');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);

figure
plot(t,emg1,'b','linewidth',1.3); hold;
plot(t,y_est1,'r','linewidth',1.5); hold;
legend('flexão canal1','FK1');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
figure
subplot(211)
plot(t,emg1,'b'); hold;
plot(t,y_est1,'r'); hold;
legend('flexão canal1','FK1');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
subplot(212)
plot(t,emg2,'b'); hold;
plot(t,y_est2,'r');hold
legend('flexão canal2','FK2');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');

figure
subplot(211)
plot(t,emg3,'b'); hold
plot(t,y_est3,'r');hold
legend('Extensão canal3','FK3');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
subplot(212)
plot(t,emg4,'b'); hold
plot(t,y_est4,'r');hold
legend('Extensão canal4','FK4');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');


figure
subplot(221)
plot(t,emg1,'b'); hold;
plot(t,y_est1,'r'); hold;
legend('flexão canal1','FK1');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);
subplot(222)
plot(t,emg2,'b'); hold;
plot(t,y_est2,'r');hold
legend('flexão canal2','FK2');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);
subplot(223)
plot(t,emg3,'b'); hold
plot(t,y_est3,'r');hold
legend('Extensão canal3','FK3');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);
subplot(224)
plot(t,emg4,'b'); hold
plot(t,y_est4,'r');hold
legend('Extensão canal4','FK4');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');ylim([0 6]);


figure
subplot(211)
plot(t,emg1,'b',t,emg2,'r'); hold;
legend('flexão canal1','flexão canal2');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
subplot(212)
plot(t,emg4,'b',t,y_est4,'r');
legend('Extensão canal3','Extensão canal4','Fusão sensorial' );
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');

figure
subplot(211)
% plot(t,emg1,'b',t,emg2,'g'); hold;
plot(t,y_est1,'b',t,y_est2,'g','linewidth',1.5);hold;
plot(t,y1,'r','linewidth',2);
legend('FK1','FK2','Fusão sensorial' );
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
subplot(212)
% plot(t,emg4,'b',t,y_est4,'g'); hold
plot(t,y_est4,'b',t,y_est3,'g','linewidth',1.5); hold ;
plot(t,y2,'r','linewidth',2);
legend('FK3','FK4','Fusão sensorial' );
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');

figure
subplot(211)
% plot(t,emg1,'b',t,emg2,'g'); hold;
plot(t,emg1,'b',t,emg2,'g','linewidth',1.5);hold;
plot(t,y1,'r','linewidth',2);
legend('canal1','canal2','Fusão sensorial' );
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
subplot(212)
% plot(t,emg4,'b',t,y_est4,'g'); hold
plot(t,emg3,'b',t,emg4,'g','linewidth',1.5); hold ;
plot(t,y2,'r','linewidth',2);
legend('canal3','canal4','Fusão sensorial' );
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');
%%
d1=[1;1];
y1=Norm(y1,5,0);
y2=Norm(y2,5,0);     

% g1=evalfis([y2' y1'],a);

 figure
 subplot(211)
 plot(t,y1,'b','linewidth',1.3); hold
 plot(t,y2,'g','linewidth',1.3);
legend('FDK flexão','FDK extensão');
ylabel('Amplitude Normalizada');
xlabel('Tempo (s)');

subplot(212)
yr=y1-y2;
yr=Norm(yr,5,0);
 plot(t,yr,'r','linewidth',1.3);
 legend('movimento Estimado');
 ylabel('Amplitude Normalizada');
 xlabel('Tempo (s)');
 %% Controle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ruído
variancia = 0.005;
x = rand(nit,1); x = 2*x-1; % Variável auxiliar
xi = sqrt(variancia*5)*x; % Ruído Branco
xi=xi;
 % modelo discretizado 
% 
a0=1; a1 =- 1.466; a2 =0.4662; 
b0 =0.006341; b1 =0.01433;
c0=1; c1=-0.2407; c2=0.0151;
% 
% a0=1; a1 =-1.0690; a2 =0.0692; 
% b0 =0.00972; b1 =0.03366;
% c0=1; c1=-0.3064; c2=0.3550;
Kp=1;
lam=(1/Kp);
A=[1 a1 a2]; B=[b0 b1 0]; C=[1 c1 c2];
TS=0.1;
GZ=tf(B,A,TS)
% Ts=Win*Ts; % tempo de amostragem
%%
%(3) DIGITAL PID CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PID tuning by "pole cancellation" and closed-loop time constant selection
    tau_mf =0.25; % given in seconds 0.25
    zd = exp(-TS/tau_mf); % desired closed-loop real discrete pole.

    s0= (1-zd)/(b0+b1);
    s1= s0*a1;
    s2= s0*a2;
    
% Sequencia de referencia
% yr=y2-y1;
% yr=Norm(yr,5,0);
% yr(1:10)=0;
nit=length(yr);

i=0; j=0;
for k=1:nit,
    % Simulated test
    if j>=((TS/Ts))
    i=i+1;
    r(i)=yr(k);        
        j=0;
    end
    j=j+1;
end
nit= length(r);
r(1:10)=0;
for k=1:2,
    ys(k)=0; us(k)=0; es(k)=0;
end
for k=3:nit,
    % Simulated test
    ys(k)= -a1*ys(k-1) -a2*ys(k-2) + b0*us(k-1) + b1*us(k-2)...
           +c0*xi(k)+c1*xi(k-2)+c2*xi(k-2);
    es(k)= r(k)-ys(k);
    us(k)=us(k-1) +(s0*es(k) +s1*es(k-1) +s2*es(k-2))/lam;
end
 tu=0:TS:nit*TS-TS; % Time vector based on Ts.

 % Ploting results
figure; % opens a new figure window
% Variancias
var_y_PID = var(ys);
var_u_PID = var(us);
ISE_PID=es*es'
% t=0:Ts:nit*Ts-Ts;

subplot(211),
    plot(tu,r,':k','linewidth',1.5); hold;
    plot(tu,ys,'r');
    title(['variância PID \sigma_\phi^2= ' num2str(var_y_PID)]);
     xlabel('tempo(s)');
     ylabel('y'); grid;
     legend('Yr','Y PID');
subplot(212),
    plot(tu,us,'r'); hold;
    title(['variância PID  \sigma_\phi^2= ' num2str(var_u_PID)]);
     xlabel('tempo(s)');
     ylabel('u'); grid;
      legend('U_ PID'); 
      
      % Controle PID digital com base no analógico paralelo
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kd= s2*Ts
Kp= -s1 -2*(Kd/Ts)
Ki= (-Kp +s0 -Kd/Ts)/Ts
td=Kd/Kp;
ti=Kp/Ki;
kc=Kp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% ARMAX
clear us y e0 ys du  
d=4;

 nit =nit-d;
p0=s0; p1=s1; p2=s2;

%o vetor H é o polinômio P vezes o C = PC 
H = conv([p0 p1 p2],C);
h0 = H(1); h1 = H(2); h2 = H(3); h3 = H(4); h4 = H(5);
%% d=4
e0=p0;
e1=p1+c1*p0+(1-a1)*e0;
e2=p2+c1*p1+c2*p0+(1-a1)*e1+(a1-a2)*e0;
e3=c1*p2+c2*p1+(1-a1)*e2+(a1-a2)*e1 +a2*e0;
f0=c2*p2+(1-a1)*e3+(a1-a2)*e2+a2*e1;
f1=(a1-a2)*e3+a2*e2
f2=a2*e3
e4=0;
e5=0;

% %  d=6
% e0=p0;
% e1=p1+c1*p0+(1-a1)*e0;
% e2=p2+c1*p1+c2*p0+(1-a1)*e1+(a1-a2)*e0;
% e3=c1*p2+c2*p1+(1-a1)*e2+(a1-a2)*e1 +a2*e0;
% e4=p2*c2+(1-a1)*e3+(a1-a2)*e2+a2*e1;
% e5=(1-a1)*e4+(a1-a2)*e3+a2*e2;
% f0=(1-a1)*e5+(a1-a2)*e4+a2*e3;
% f1=(a1-a2)*e5+a2*e4;
% f2=a2*e5;

Q= lam*[1 0 0 0 0];
BE=[conv(B(1:2),[e0 e1 e2 e3 e4 e5])]
CQ=conv(C,Q)
BECQ=BE+CQ;
Ax = conv([1 -1],A);
p=[p0 p1 p2];
BP=[conv(B(1:2),p)];
AxQ=[conv(Ax,Q)];
r(1:10)=0;
for k=1:7,
    u(k)=0; y(k)=0; e(k)=0;  du(k)=0;
end 

for k=7:nit,
    y(k)= -a1*y(k-1) -a2*y(k-2) + b0*u(k-1) + b1*u(k-2)...
           +c0*xi(k)+c1*xi(k-2)+c2*xi(k-2);%+... 
%            +v(k) +a1*v(k-1) +a2*v(k-2); %pertubação

       
%     // GMV
    du(k) = (1/BECQ(1))*( -BECQ(2)*du(k-1)-BECQ(3)*du(k-2)-BECQ(4)*du(k-3)-BECQ(5)*du(k-4)-BECQ(6)*du(k-5)-BECQ(7)*du(k-6) ...
             +h0*r(k+d) +h1*r(k+d-1) +h2*r(k+d-2) +h3*r(k+d-3) +h4*r(k+d-4) ...
             -f0*y(k) -f1*y(k-1) -f2*y(k-2) );
    
      e=r(k)-y(k);
      u(k)=u(k-1) +du(k);
%   u(k)=du(k);
end 
ISE_armax=e*e'
% Plotando os resultados
% t=0:Ts:nit*Ts-Ts;
%scf(); t=0:Ts:nit*Ts-Ts;
 figure;
% Variancias
var_y_armax = var(y);
var_u_armax = var(u);

% t=0:Ts:nit*Ts-Ts;
tu=0:TS:nit*TS-TS;
subplot(211),
    plot(tu,r(1:nit),':k','linewidth',1.5); hold;
    plot(tu,y,'r');
    y_mx=y;
    title(['variância GMV_ ARMAX  \sigma_\phi^2= ' num2str(var_y_armax)]);
     xlabel('tempo(s)');
     ylabel('y'); grid;
     legend('Yr','Y GMV_ ARMAX');
subplot(212),
    u_mx=u;
    plot(tu,u,'r'); hold;
    title(['variância GMV_.ARMAX \sigma_\phi^2= ' num2str(var_u_armax)]);
     xlabel('tempo(s)');
     ylabel('u'); grid;
      legend('U_ ARMAX'); 
%% Análise GMV
z = tf([0 0 0 0 1],1,TS,'variable','z^-1');
% z=1;
BP=tf(BP,1,TS,'variable','z^-1');
CQ=tf(CQ,1,TS,'variable','z^-1');
BE=tf(BE,1,TS,'variable','z^-1');
AxQ=tf(AxQ,1,TS,'variable','z^-1');

S_gmv=(BE+CQ)/(BP+AxQ);
SC_gmv=(BP)/(BP+AxQ);

S_pid=(BP*z)/(AxQ+BP*z);
SC_pid=(CQ)/(AxQ+BP*z);

ganhos=dsigma(SC_gmv,TS);
Mt = max(ganhos); clear ganhos;
ganhos=dsigma(S_gmv,TS);
Ms = max(ganhos); clear ganhos;
         %MG   
MG_gmv = max( (Ms/(Ms-1)), (1+ (1/Mt)) );
MGgmv_dB = 20*log10(MG_gmv)
        % MF
MF_gmv = max( (2*asin(1/(2*Ms))*(180/pi)), ...
               (2*asin(1/(2*Mt))*(180/pi)) )

figure
        dsigma(SC_gmv,TS)
        hold
        dsigma(S_gmv,TS)
%         title(['{GMV: MG >= }' num2str(MGgmv_dB) 'dB' ...
%           '{ e MF >= }' num2str(MF_gmv) '{\circ}.']);
          legend('|T_GMV|','|S_GMV|');
         ylabel('Ganhos Principais');
         xlabel('Frequencia (Hz)');
% Analise PID
% figure
% S_pid=(BP)/(AxQ+BP);
% SC_pid=(CQ)/(AxQ+BP);
S_pid=(BP*z)/(AxQ+BP*z);
SC_pid=(CQ)/(AxQ+BP*z);
[ganhosSC wSC]=dsigma(SC_pid,TS);
Mt = max(ganhosSC); 
[ganhosS wS]=dsigma(S_pid,TS);
Ms = max(ganhosS); 
         %MG   
MG_pid = max( (Ms/(Ms-1)), (1+ (1/Mt)) );
MGpid_dB = 20*log10(MG_pid)
        % MF
MF_pid = max( (2*asin(1/(2*Ms))*(180/pi)), ...
               (2*asin(1/(2*Mt))*(180/pi)) )

        dsigma(SC_pid,TS)
%         semilogx(wSC,20*log10(ganhosSC),'k','linewidth',2)
%         hold
        dsigma(S_pid,TS)
%         semilogx(wS,20*log10(ganhosS),'k','linewidth',2)
        
%         title(['{pid: MG >= }' num2str(MGpid_dB) 'dB' ...
%           '{ e MF >= }' num2str(MF_pid) '{\circ}.']);
        title(' ')
          legend('|GMV Sens-Co|','|GMV Sens|','|PID Sens-Co|','|PID Sens|');
         ylabel('Ganhos Principais');
         xlabel('Frequencia');
         
         set(gca,'fontsize',12);
         set(gca,'linewidth',1);