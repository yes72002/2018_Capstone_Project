clear
close all
clc
%{
Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 1000;             % Length of signal
t = (0:L-1)*T;        % Time vector
S = 0.3*sin(2*pi*50*t) ;
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%}
time=0.1;
fsw=3000;
f=60;
n=fsw/f;
N=time*fsw;
freqStep = fsw/N;
freq = freqStep*(0:N/2-1);
t = 0 : 1/fsw :time; 
r = sin(2*pi*f*t); 
fft_r = fft(r,N);
mag_r = abs(fft_r)/(N/2); 
a=mag_r(1:1:N/2);
plot(freq,a)
