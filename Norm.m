function y_n=Norm(x,amplitude_MAX, amplitude_MIN);
%% Normaliza��o de dados
M=max(x);
M1=min(x);
MM=max([M1;M]);
x=(x -M1)*(amplitude_MAX-amplitude_MIN)/(M-M1) + amplitude_MIN;
y_n=x;
