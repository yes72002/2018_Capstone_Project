clear%6/29利用w=0.5產生方波
close all
clc
time = 0.1;
fsw = 5000;
bit = 100;
amp = 0.5;
f = 50;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t); 
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3); 
% subplot(4,1,1),plot(t,s1,t,s2,t,s3);

w = 0.5;
duty1=s1+w;
out1=zeros(1,fsw*bit*time);
for i=1:1:500
x=bit*(i-1)+1:1:bit*(i-1)+fix((bit-bit*duty1(i))/2);
out1(x)=0;
x=bit*(i-1)+fix((bit-bit*duty1(i))/2+1):1:bit*(i-1)+fix((bit+bit*duty1(i))/2);
out1(x)=1;
x=bit*(i-1)+fix((bit+bit*duty1(i))/2+1):1:bit*(i-1)+bit;
out1(x)=0;
end
duty2=s2+w;
out2=zeros(1,fsw*bit*time);
for i=1:1:500
x=bit*(i-1)+1:1:bit*(i-1)+fix((bit-bit*duty2(i))/2);
out2(x)=0;
x=bit*(i-1)+fix((bit-bit*duty2(i))/2+1):1:bit*(i-1)+fix((bit+bit*duty2(i))/2);
out2(x)=1;
x=bit*(i-1)+fix((bit+bit*duty2(i))/2+1):1:bit*(i-1)+bit;
out2(x)=0;
end
duty3=s3+w;
out3=zeros(1,fsw*bit*time);
for i=1:1:500
x=bit*(i-1)+1:1:bit*(i-1)+fix((bit-bit*duty3(i))/2);
out3(x)=0;
x=bit*(i-1)+fix((bit-bit*duty3(i))/2+1):1:bit*(i-1)+fix((bit+bit*duty3(i))/2);
out3(x)=1;
x=bit*(i-1)+fix((bit+bit*duty3(i))/2+1):1:bit*(i-1)+bit;
out3(x)=0;
end

% subplot(4,1,2),plot(out1);
% subplot(4,1,3),plot(out2);
% subplot(4,1,4),plot(out3);
%{
N=time*fsw*bit;
freqStep = 10;
freq = freqStep*(0:N/2-1);
fft_out1 = fft(out1,N);
mag_out1 = abs(fft_out1)/(length(fft_out1)/2); 
a1=mag_out1(1:1:N/2);
fft_out2 = fft(out2,N);
mag_out2 = abs(fft_out2)/(length(fft_out2)/2); 
a2=mag_out2(1:1:N/2);
fft_out3 = fft(out3,N);
mag_out3 = abs(fft_out3)/(length(fft_out3)/2); 
a3=mag_out3(1:1:N/2);
for i=2:fsw*bit*time/2
  freq(i-1)=freq(i);
  a1(i-1)=a1(i);
  a2(i-1)=a2(i);
  a3(i-1)=a3(i);
end
figure;
subplot(3,1,1),plot(freq,a1);
subplot(3,1,2),plot(freq,a2);
subplot(3,1,3),plot(freq,a3)

harm1 =a1( 12: 6: end/2 ); 
THD1 = 100*sqrt(  sum( harm1.^2) /a1(6)^2 );
harm2 = a2( 12: 6: end/2 ); 
THD2 = 100*sqrt(  sum( harm2.^2) /a2(6)^2 );
harm3 = a3( 12: 6: end/2 ); 
THD3 = 100*sqrt(  sum( harm3.^2) /a3(6)^2 );
%}

out=[out1;out2;out3];
N = time*fsw*bit;
Van=[ 2/3 -1/3 -1/3]*out;
Vbn=[-1/3  2/3 -1/3]*out;
Vcn=[-1/3 -1/3  2/3]*out;
fft_Van=fft(Van,N)/N;%6:0.499且兩端都有值
mag_Van=abs(fft_Van)*2;
fft_Vbn=fft(Vbn,N)/N;
mag_Vbn=abs(fft_Vbn)*2;
fft_Vcn=fft(Vcn,N)/N;
mag_Vcn=abs(fft_Vcn)*2;
i=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:到中間
mag_Van_1=mag_Van(1,i);
mag_Vbn_1=mag_Vbn(1,i);
mag_Vcn_1=mag_Vcn(1,i);
% figure;%=================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(倍頻^2/基頻^2)
THD_Van=100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n',THD_Vcn);
%ITHD=頻率除以本身值
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
for i=1:fsw*bit*0.02/2-1
    mag_ithd_out1(i)=mag_Van_1(i)/(10*i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(10*i); 
    mag_ithd_out3(i)=mag_Vcn_1(i)/(10*i);  
end
ITHD_Van=100*((sum(mag_ithd_out1.^2)-(mag_ithd_out1(1,1).^2))/(mag_ithd_out1(1,1).^2)).^(1/2);
ITHD_Vbn=100*((sum(mag_ithd_out2.^2)-(mag_ithd_out2(1,1).^2))/(mag_ithd_out2(1,1).^2)).^(1/2);
ITHD_Vcn=100*((sum(mag_ithd_out3.^2)-(mag_ithd_out3(1,1).^2))/(mag_ithd_out3(1,1).^2)).^(1/2);
fprintf('ITHD1=%f\n',ITHD_Van);
fprintf('ITHD2=%f\n',ITHD_Vbn);
fprintf('ITHD3=%f\n',ITHD_Vcn);
%該算的,此行以上都算出來了,接下來是畫圖
n2=2:1:N/2;%去頭(頭第一個值有誤為零)
mag_Van_2(1,n2-1)=mag_Van(1,n2);%5:0.499,10:0.00007,15...
mag_Vbn_2(1,n2-1)=mag_Vbn(1,n2); 
mag_Vcn_2(1,n2-1)=mag_Vcn(1,n2); 
mag_Van_3=zeros(1,(N/2-1)*10);
mag_Vbn_3=zeros(1,(N/2-1)*10);
mag_Vcn_3=zeros(1,(N/2-1)*10);
for k=1:1:N/2-1;%5->50
mag_Van_3(1,(bit*fsw/N)*k)=mag_Van_2(1,k);%50:0.499
mag_Vbn_3(1,(bit*fsw/N)*k)=mag_Vbn_2(1,k);
mag_Vcn_3(1,(bit*fsw/N)*k)=mag_Vcn_2(1,k);
end
% figure;
% semilogx(mag_Van_3,'r'),title('mag Van 3');
% end_time=clock;
% execution_time=end_time-start_time;