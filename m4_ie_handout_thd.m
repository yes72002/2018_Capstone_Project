clear
close all
clc
t = 0: 1/5e3 : 0.1-1/5e3; 
r_temp = [ones( 1, 5e3/100 ) -ones( 1, 5e3/100 )];
r = [r_temp r_temp r_temp r_temp r_temp];
figure; 
plot(t,r)

fft_r = fft(r,length(r)); 
mag_r = abs(fft_r)/(length(fft_r)/2); 
freq = (0:length(fft_r)/2-1)*5e3/length(fft_r); 
%freq = (0:249)*10; 
figure; 
plot(freq,mag_r(1:end/2));
figure;
semilogx(mag_r(1:end/2));

harm = mag_r( 11: 5: end/2 ); 
THD = 100*sqrt(  sum( harm.^2) /mag_r(6)^2 );