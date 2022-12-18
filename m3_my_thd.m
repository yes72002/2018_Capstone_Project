clear
close all
clc
t = 0: 1/5e3 : 0.1 - 1/5e3;
r_temp = [ ones(1,5e3/100) -ones(1,5e3/100) ]; 
r = [r_temp r_temp r_temp r_temp r_temp];
figure;
plot( t, r )

N = 500;
s = fft(r,N);      
freqStep = 5000/N;
mag_r = abs(s)/(N/2); 
freq = freqStep*(0:N/2-1);
a = mag_r(1:1:N/2);
figure;
plot(freq,a)

%THD=100*sqar{(I2^2+I3^2+I4^2+...)/I1^2}
%s(6)--s(496)246
up=0;
for i=1:24
    up=real(a(10*i+6)^2)+up;  %a=s(16)^2+s(26)^2+s(36)^2;
end
down=real(a(6)^2);
THD=100*sqrt(up/down);



