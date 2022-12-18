clear%6/30利用w產生方波
close all
clc
time = 0.1;
fsw = 5000;
bit = 300;
amp = 0.5;
f = 50;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t); 
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3); 
figure;
subplot(5,1,1),plot(t,s1,t,s2,t,s3);
ALLS=[s1;s2;s3];
sinmax=max(ALLS);
sinmin=min(ALLS);
ratio1=1;
ratio0=1;
w1=1-sinmax;
w0=-sinmin;
w=w1*ratio1/(ratio1+ratio0)+w0*ratio0/(ratio1+ratio0);
duty1=s1+w;
out1=zeros(1,fsw*bit*time);
for i=1:1:fsw*time
x=bit*(i-1)+1:1:bit*(i-1)+fix((bit-bit*duty1(i))/2);
out1(x)=0;
x=bit*(i-1)+fix((bit-bit*duty1(i))/2+1):1:bit*(i-1)+fix((bit+bit*duty1(i))/2);
out1(x)=1;
x=bit*(i-1)+fix((bit+bit*duty1(i))/2+1):1:bit*(i-1)+bit;
out1(x)=0;
end
duty2=s2+w;
out2=zeros(1,fsw*bit*time);
for i=1:1:fsw*time
x=bit*(i-1)+1:1:bit*(i-1)+fix((bit-bit*duty2(i))/2);
out2(x)=0;
x=bit*(i-1)+fix((bit-bit*duty2(i))/2+1):1:bit*(i-1)+fix((bit+bit*duty2(i))/2);
out2(x)=1;
x=bit*(i-1)+fix((bit+bit*duty2(i))/2+1):1:bit*(i-1)+bit;
out2(x)=0;
end
duty3=s3+w;
out3=zeros(1,fsw*bit*time);
for i=1:1:fsw*time
x=bit*(i-1)+1:1:bit*(i-1)+fix((bit-bit*duty3(i))/2);
out3(x)=0;
x=bit*(i-1)+fix((bit-bit*duty3(i))/2+1):1:bit*(i-1)+fix((bit+bit*duty3(i))/2);
out3(x)=1;
x=bit*(i-1)+fix((bit+bit*duty3(i))/2+1):1:bit*(i-1)+bit;
out3(x)=0;
end
out=[out1;out2;out3];
subplot(5,1,2),plot(out1);
subplot(5,1,3),plot(out2);
subplot(5,1,4),plot(out3);

plusw1=out1;%給Untitled12比較用
plusw2=out2;
plusw3=out3;

% out4=zeros(1,fsw*bit*time);
% for i=1:1:fsw*bit*time%產生111,110,100,000的方波(震幅為2,1,0)
%     if out1(i)==0 && out2(i)==0 && out3(i)==0
%         out4(i)=0;
%     elseif out1(i)==1 && out2(i)==1 && out3(i)==1
%         out4(i)=2;
%     else
%         out4(i)=1;
%     end
% end
% subplot(5,1,5),plot(out4);

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
% subplot(3,1,1),plot(freq,a1)
% subplot(3,1,2),plot(freq,a2);
% subplot(3,1,3),plot(freq,a3)
plot(freq,a1);
figure;
semilogx(a1);

k=f/freqStep;
harm1 = a1(2*k:k:end/2); 
THD1 = 100*sqrt(sum( harm1.^2) /a1(k)^2);
harm2 = a2(2*k:k:end/2); 
THD2 = 100*sqrt(sum( harm2.^2) /a2(k)^2);
harm3 = a3(2*k:k:end/2); 
THD3 = 100*sqrt(sum( harm3.^2) /a3(k)^2);

A = [ 2/3 -1/3 -1/3 ; -1/3 2/3 -1/3 ; -1/3 -1/3 2/3];  
B = [ out1 ; out2 ; out3 ];  
C = 48*A*B  ; 
figure;
plot(C(1,1:1:fsw*bit*time));
hold on;
plot(C(2,1:1:fsw*bit*time));
hold on;
plot(C(3,1:1:fsw*bit*time));
hold on;
% figure
% semilogx(C(1,1:1:fsw*bit*time))