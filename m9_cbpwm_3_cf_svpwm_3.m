clear%7/7¤ñ¸û
close all
clc
run('m7_cbpwm_changed.m');

time = 0.1;
fsw = 5000;
bit = 300;
amp = 0.5;
f = 60;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*sin(2*pi*f*t); 
s2 = amp*sin(2*pi*f*t+2*pi/3); 
s3 = amp*sin(2*pi*f*t+4*pi/3); 
subplot(4,1,1),plot(t,s1,t,s2,t,s3);

mapping=2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1 1 1];
sindq=mapping*[s1;s2;s3];
sindqxy=[sindq(1,:);sindq(2,:)];
sindqx=sindqxy(1,:);
sindqy=sindqxy(2,:);
i=1:1:500;
angle(i)=atan2(sindqy(i),sindqx(i));
M4=[2/3 ; 0];
M6=[1/3 ; sqrt(3)/3];
M2=[-1/3 ; sqrt(3)/3];
M3=[-2/3 ; 0];
M1=[-1/3 ; -sqrt(3)/3];
M5=[1/3 ; -sqrt(3)/3];
duty=zeros(2,fsw*time);
sector=zeros(1,500);
a0=zeros(1,500);
a7=zeros(1,500);
aless=zeros(1,500);
amore=zeros(1,500);
out1=zeros(1,fsw*bit*time);
out2=zeros(1,fsw*bit*time);
out3=zeros(1,fsw*bit*time);
for i=1:fsw*time
if angle(i)>0 && angle(i)<=pi/3;%sector1
    sector(i)=1;
    duty=[M4 M6]\[sindqx(i) ; sindqy(i)];
    a4=duty(1);
    a6=duty(2);
    a0(i)=(1-a4-a6)/2;
    a7(i)=a0(i);
    aless(i)=a4;
    amore(i)=a6;
end
if angle(i)>-pi/3 && angle(i)<=0;%sector6
    sector(i)=6;
    duty=[M4 M5]\[sindqx(i) ; sindqy(i)];
    a4=duty(1);
    a5=duty(2);
    a0(i)=(1-a4-a5)/2;
    a7(i)=a0(i);
    aless(i)=a4;
    amore(i)=a5;
end
if angle(i)>pi/3 && angle(i)<=2*pi/3;%sector2
    sector(i)=2;
    duty=[M2 M6]\[sindqx(i) ; sindqy(i)];
    a2=duty(1);
    a6=duty(2);
    a0(i)=(1-a2-a6)/2;
    a7(i)=a0(i);
    aless(i)=a2;
    amore(i)=a6;
end   
if angle(i)>-2*pi/3 && angle(i)<=-pi/3;%sector5
    sector(i)=5;
    duty=[M1 M5]\[sindqx(i) ; sindqy(i)];
    a1=duty(1);
    a5=duty(2);
    a0(i)=(1-a1-a5)/2;
    a7(i)=a0(i);
    aless(i)=a1;
    amore(i)=a5;
end  
if angle(i)>2*pi/3 && angle(i)<=pi;%sector3
    sector(i)=3;
    duty=[M2 M3]\[sindqx(i) ; sindqy(i)];
    a2=duty(1);
    a3=duty(2);
    a0(i)=(1-a2-a3)/2;
    a7(i)=a0(i);
    aless(i)=a2;
    amore(i)=a3;
end
if angle(i)>-pi && angle(i)<=-2*pi/3;%sector4
    sector(i)=4;
    duty=[M1 M3]\[sindqx(i) ; sindqy(i)];
    a1=duty(1);
    a3=duty(2);
    a0(i)=(1-a1-a3)/2;
    a7(i)=a0(i);
    aless(i)=a1;
    amore(i)=a3;
end
end
for i=1:fsw*time
    switch sector(i)
        case 1
            j=bit*(i-1)+1:1:bit*(i-1)+fix(a0(i)*bit/2);
            out1(j)=0;
            j=bit*(i-1)+fix(a0(i)*bit/2)+1:1:bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit);
            out1(j)=1;
            j=bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out1(j)=0;
    
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i))*bit/2);
            out2(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit);
            out2(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out2(j)=0;
    
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2);
            out3(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit);
            out3(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out3(j)=0;

        case 2
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i))*bit/2);
            out1(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit);
            out1(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out1(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix(a0(i)*bit/2);
            out2(j)=0;
            j=bit*(i-1)+fix(a0(i)*bit/2)+1:1:bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit);
            out2(j)=1;
            j=bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out2(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2);
            out3(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit);
            out3(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out3(j)=0;
        case 3
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2);
            out1(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit);
            out1(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out1(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix(a0(i)*bit/2);
            out2(j)=0;
            j=bit*(i-1)+fix(a0(i)*bit/2)+1:1:bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit);
            out2(j)=1;
            j=bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out2(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i))*bit/2);
            out3(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit);
            out3(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out3(j)=0;
        case 4
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2);
            out1(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit);
            out1(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out1(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i))*bit/2);
            out2(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit);
            out2(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out2(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix(a0(i)*bit/2);
            out3(j)=0;
            j=bit*(i-1)+fix(a0(i)*bit/2)+1:1:bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit);
            out3(j)=1;
            j=bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out3(j)=0;
        case 5 
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i))*bit/2);
            out1(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit);
            out1(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out1(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2);
            out2(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit);
            out2(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out2(j)=0;
            
            j=bit*(i-1)+1:1:bit*(i-1)+fix(a0(i)*bit/2);
            out3(j)=0;
            j=bit*(i-1)+fix(a0(i)*bit/2)+1:1:bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit);
            out3(j)=1;
            j=bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out3(j)=0;
        case 6  
            j=bit*(i-1)+1:1:bit*(i-1)+fix(a0(i)*bit/2);
            out1(j)=0;
            j=bit*(i-1)+fix(a0(i)*bit/2)+1:1:bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit);
            out1(j)=1;
            j=bit*(i-1)+fix((a0(i)/2+aless(i)+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out1(j)=0;
    
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2);
            out2(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i)+amore(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit);
            out2(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i)+amore(i))/2+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out2(j)=0;
    
            j=bit*(i-1)+1:1:bit*(i-1)+fix((a0(i)+aless(i))*bit/2);
            out3(j)=0;
            j=bit*(i-1)+fix((a0(i)+aless(i))*bit/2)+1:1:bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit);
            out3(j)=1;
            j=bit*(i-1)+fix(((a0(i)+aless(i))/2+amore(i)+a7(i))*bit)+1:1:bit*(i-1)+bit;
            out3(j)=0;
    end
end
subplot(4,1,2),plot(out1);
subplot(4,1,3),plot(out2);
subplot(4,1,4),plot(out3);

sususu1=out1;
sususu2=out2;
sususu3=out3;
for i=1:1:fsw*time*bit
   cf1(i)=plusw1(i)-sususu1(i);
   cf2(i)=plusw2(i)-sususu2(i);  
   cf3(i)=plusw3(i)-sususu3(i);
end
diff=sum(abs(cf1))+sum(abs(cf2))+sum(abs(cf3))
subplot(3,1,1),plot(cf1);
subplot(3,1,2),plot(cf2);
subplot(3,1,3),plot(cf3);

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

harm1 = a1( 12: 6: end/2 ); 
THD1 = 100*sqrt(  sum( harm1.^2) /a1(6)^2 );
harm2 = a2( 12: 6: end/2 ); 
THD2 = 100*sqrt(  sum( harm2.^2) /a2(6)^2 );
harm3 = a3( 12: 6: end/2 ); 
THD3 = 100*sqrt(  sum( harm3.^2) /a3(6)^2 );

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
