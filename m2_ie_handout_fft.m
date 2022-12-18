clear
close all
clc
time = 0 : 1/3000 : 0.1-1/3000; 
r = sin(2*pi*60*time); 
fft_r = fft(r,300);
mag_r = abs(fft_r)/(length(fft_r)/2); 
f=0:10:3000-10;
plot(f,mag_r)
