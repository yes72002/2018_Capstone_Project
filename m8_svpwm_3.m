clear%7/4and7/5利用duty產生方波
close all
clc
start_time=clock;
time = 0.1;
fsw = 5000;
bit = 400;
amp = 0.5;
f = 50;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t); 
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3); 
a = 2*pi/3;
mapping=2/3*[cos(0) cos(a) cos(2*a);sin(0) sin(a) sin(2*a);1 1 1];
sindq=mapping*[s1;s2;s3];
sindqxy=[sindq(1,:);sindq(2,:)];
sindqx=sindqxy(1,:);
sindqy=sindqxy(2,:);
angle=atan2(sindqy,sindqx);
for i=1:1:fsw*time
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
%plot(sindqy,sindqx);
M4=[2/3 ; 0];
M6=[1/3 ; sqrt(3)/3];
M2=[-1/3 ; sqrt(3)/3];
M3=[-2/3 ; 0];
M1=[-1/3 ; -sqrt(3)/3];
M5=[1/3 ; -sqrt(3)/3];

duty=zeros(2,fsw*time);
sector=zeros(1,fsw*time);
aless=zeros(1,fsw*time);
amore=zeros(1,fsw*time);
v1=zeros(3,fsw*time);
v2=zeros(3,fsw*time);
for i=1:fsw*time%判斷sector
if angle(i)>0 && angle(i)<=pi/3;%sector1
    sector(i)=1;
    duty=[M4 M6]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    v1(:,i)=[1;0;0];
    v2(:,i)=[1;1;0];
end
if angle(i)>pi/3 && angle(i)<=2*pi/3;%sector2
    sector(i)=2;
    duty=[M2 M6]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    v1(:,i)=[0;1;0];
    v2(:,i)=[1;1;0];
end   
if angle(i)>2*pi/3 && angle(i)<=pi;%sector3
    sector(i)=3;
    duty=[M2 M3]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    v1(:,i)=[0;1;0];
    v2(:,i)=[0;1;1];
end
if angle(i)>pi && angle(i)<=4*pi/3;%sector4
    sector(i)=4;
    duty=[M1 M3]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    v1(:,i)=[0;0;1];
    v2(:,i)=[0;1;1];
end
if angle(i)>4*pi/3 && angle(i)<=5*pi/3;%sector5
    sector(i)=5;
    duty=[M1 M5]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    v1(:,i)=[0;0;1];
    v2(:,i)=[1;0;1];
end 
if angle(i)>5*pi/3 && angle(i)<=2*pi;%sector6
    sector(i)=6;
    duty=[M4 M5]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    v1(:,i)=[1;0;0];
    v2(:,i)=[1;0;1];
end
end
az=1-aless-amore;
a0=az/2;
a7=az/2;

T1=aless;%無限跟有限精準度
T2=amore;
Tz=az;
T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
Tzbit=bit/2-T1bit-T2bit;
TTT=T1bit+T2bit+Tzbit;
%{
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
%}
for i=1:fsw*time%主要的大程式碼
    out000=repmat([0;0;0],1,floor(Tzbit(i)/2));
    out001=repmat(v1(:,i),1,T1bit(i));
    out011=repmat(v2(:,i),1,T2bit(i));
    out111=repmat([1;1;1],1,ceil(Tzbit(i)/2));
    out0127=[out000 out001 out011 out111 out111 out011 out001 out000];
    out(:,bit*(i-1)+1:bit*i)=out0127;
end
out2=circshift(out,[0 1]);
out3=[out;out2];
diff1=sum(abs(out(1,:)-out2(1,:)));
diff2=sum(abs(out(2,:)-out2(2,:)));
diff3=sum(abs(out(3,:)-out2(3,:)));
fprintf('out1切換次數=%f\n',diff1);
fprintf('out2切換次數=%f\n',diff2);
fprintf('out3切換次數=%f\n',diff3);

figure;
subplot(4,1,1),plot(t,s1,t,s2,t,s3);
title('三相輸入');
subplot(4,1,2),plot(out(1,:));
title('out1輸出');
subplot(4,1,3),plot(out(2,:));
title('out2輸出');
subplot(4,1,4),plot(out(3,:));
title('out3輸出');

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
%}
N = time*fsw*bit;
%{
Van=[ 2/3 -1/3 -1/3]*out;
Vbn=[-1/3  2/3 -1/3]*out;
Vcn=[-1/3 -1/3  2/3]*out;
fft_out1 = fft(Van,N); 
mag_out1 = abs(fft_out1)/(length(fft_out1)/2); 
fft_out2 = fft(Vbn,N); 
mag_out2 = abs(fft_out2)/(length(fft_out2)/2); 
fft_out3 = fft(Vcn,N); 
mag_out3 = abs(fft_out3)/(length(fft_out3)/2); 
freq = (1:length(fft_out1)/2-1)*fsw*bit/length(fft_out1); 
%freq = (0:150000/2-1)*10; 
draw1=mag_out1(2:1:N/2);%為了畫fft才做這樣的轉換
draw2=mag_out2(2:1:N/2);
draw3=mag_out3(2:1:N/2);

figure;
% subplot(5,1,1),plot(freq,draw1);
% subplot(5,1,2),plot(freq,draw2);
% subplot(5,1,3),plot(freq,draw3);
plot(freq,draw1);
title('out1 fft');
figure;
semilogx(mag_out1(1:end/2));
title('semilogx out1 fft');
k=f*time;
harm1 = mag_out1( 2*k+1: k: end/2 );
THD1 = 100*sqrt(  sum( harm1.^2)/mag_out1(k+1)^2 );
harm2 = mag_out2( 2*k+1: k: end/2 );
THD2 = 100*sqrt(  sum( harm2.^2)/mag_out2(k+1)^2 );
harm3 = mag_out3( 2*k+1: k: end/2 );
THD3 = 100*sqrt(  sum( harm3.^2)/mag_out3(k+1)^2 );
mag_ithd_out1=zeros(1,fsw*bit*time/2);
mag_ithd_out2=zeros(1,fsw*bit*time/2);
mag_ithd_out3=zeros(1,fsw*bit*time/2);
for i=1:fsw*bit*time
    mag_ithd_out1(i)=mag_out1(i)/(10*i); 
    mag_ithd_out2(i)=mag_out2(i)/(10*i); 
    mag_ithd_out3(i)=mag_out3(i)/(10*i); 
end
harm_ithd1 = mag_ithd_out1( 2*k+1: k: end/2 );
ITHD1 = 100*sqrt(  sum( harm_ithd1.^2)/mag_ithd_out1(k+1)^2 );
harm_ithd2 = mag_ithd_out2( 2*k+1: k: end/2 );
ITHD2 = 100*sqrt(  sum( harm_ithd2.^2)/mag_ithd_out2(k+1)^2 );
harm_ithd3 = mag_ithd_out3( 2*k+1: k: end/2 );
ITHD3 = 100*sqrt(  sum( harm_ithd3.^2)/mag_ithd_out3(k+1)^2 );
fprintf('VTHD1=%f\n',THD1);
fprintf('VTHD2=%f\n',THD2);
fprintf('VTHD3=%f\n\n',THD3);
fprintf('ITHD1=%f\n',ITHD1);
fprintf('ITHD2=%f\n',ITHD2);
fprintf('ITHD3=%f\n\n',ITHD3);

A = [ 2/3 -1/3 -1/3 ; -1/3 2/3 -1/3 ; -1/3 -1/3 2/3];  
B = [ out(1,:) ; out(2,:) ; out(3,:) ];  
C = 48*A*B  ; 
figure;
plot(C(1,1:1:fsw*bit*time));
hold on;
plot(C(2,1:1:fsw*bit*time));
hold on;
plot(C(3,1:1:fsw*bit*time));
title('相電壓');
% figure
% semilogx(C(1,1:1:fsw*bit*time))    
end_time=clock;
execution_time=end_time-start_time;
%}
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
figure;
semilogx(mag_Van_1,'r'),title('mag Van 1');
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
figure;
semilogx(mag_Van_3,'r'),title('mag Van 3');
end_time=clock;
execution_time=end_time-start_time;
ITHDavg=(ITHD_Van+ITHD_Vbn+ITHD_Vcn)/3;