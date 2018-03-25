function analise_dos_dados(t,y1);
% clear all;
% clc; close all;
% [X] = importdata('teste3.TXT');
% [t, y1, y2, y3, y4] = textread('teste1.TXT', ...
% '%f %f %f %f %f',60001);
% 
%% Análises
Ts=5.0000e-04;
A = y1;
[N,~] = size(A);
time = A(N);
N = N-1;
X = A(1 : N);
Fs = N/time;
xdft = fft(X);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

figure
T=(0:Ts: length(y1)*Ts -Ts);
plot(T,A);
grid on
title('Sinal Eletromiográfico do Grupo Flexor durante a Flexão do punho')
xlabel('Tempo (ms)')
ylabel('Tensão (mV)')

figure
plot(freq,10*log10(psdx));
grid on
title('Periodogram Using FFT')
xlabel('Frequência (Hz)')
ylabel('Potência/Frequencia (dB/Hz)')



%% Loading Raw Signal of EMG %%

emg = y1;

%% Removing the DC Offset from the Signal %%

% emg2 = detrend(emg); % detrend removes the mean value or linear trend from a vector or matrix. 
Fs = 1/Ts; Fn = Fs/2; Wn = 100/Fn; 
[b,a]=butter(1,Wn,'HIGH'); % low pass filter to cut off frequency at 10Hz, using a 5th order filter and a samppling frequencyt of 1000Hz.
emg2=filtfilt(b,a,emg);

 %% Rectifcication of the Signal 
 
 emg3=abs(emg2);
 
%% Creating the Envelope of the EMG Signal %%

% emg2 = detrend(emg); % detrend removes the mean value or linear trend from a vector or matrix. 
Fs = 1/Ts; Fn = Fs/2; Wn = 5/Fn; 
[b,a]=butter(1,Wn,'low'); % low pass filter to cut off frequency at 10Hz, using a 5th order filter and a samppling frequencyt of 1000Hz.
emg4=filtfilt(b,a,emg3);

%% Creating the FFT %%

L = length(emg);
% Ts = mean(diff(t));                                 % Sampling Time
Fs = 1/Ts;                                          % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency

FEMG = fft(emg)*2/L;                                % Fourier Transform (Normalised)
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
Iv = 1:length(Fv);                                  % Index Vector

figure
plot(Fv, abs(FEMG(Iv)))
grid
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% Plotting all 4 Plots %% 
emg1=y1;
figure
subplot(2,2,1)       
plot(t,emg1,'k')
grid on;
xlabel ('Tempo (s)');
ylabel ('Amplitude (uV)'); 
title ('Sinal EMG Bruto')

subplot(2,2,2)      
plot (t,emg2,'k')        
grid on;
xlabel ('Tempo (s)');
ylabel ('Amplitude (uV)');
title('Sinal sem DC Offset')

subplot(2,2,3)       
plot(t,emg3,'k')
grid on; 
xlabel('Tempo (s)');
ylabel ('Amplitude (uV)'); 
title ('Sinal EMG Retificado');
 

subplot(2,2,4)       
plot(t,emg4,'k')
grid on;
xlabel ('Tempo (s)');
ylabel ('Amplitude (uV)');
title ('Envelope linear do sinal EMG');
%%
% EMG =(y1);
% L = length(t);
% Ts = mean(diff(t));                                 % Sampling Time
% % Ts=5.0000e-04;
% Fs = 1/Ts;                                          % Sampling Frequency
% Fn = Fs/2;                                          % Nyquist Frequency
% FEMG = fft(EMG)*2/L;                                % Fourier Transform (Normalised)
% Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
% Iv = 1:length(Fv);                                  % Index Vector
% figure
% plot(Fv, abs(FEMG(Iv)))
% grid
end