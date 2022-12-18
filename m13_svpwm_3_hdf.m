clear%7/20,7/24利用最低的HDF寫出序列
close all
clc
start_time=clock;
time = 0.02;
fsw = 15000;
bit = 400;
amp = 0.5;
f = 50;
n = fsw*time;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t);
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3); 
a = 2*pi/3;
mapping=2/3*[cos(0) cos(a) cos(2*a);sin(0) sin(a) sin(2*a);1 1 1];
sindqxy=mapping*[s1;s2;s3];
sindqx=sindqxy(1,:);
sindqy=sindqxy(2,:);
angle=atan2(sindqy,sindqx);
for i=1:1:fsw*time
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
A = [2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3];
DQ1=mapping*A*[0;0;1];
DQ2=mapping*A*[0;1;0];
DQ3=mapping*A*[0;1;1];
DQ4=mapping*A*[1;0;0];
DQ5=mapping*A*[1;0;1];
DQ6=mapping*A*[1;1;0];
M1=[DQ1(1);DQ1(2)];
M2=[DQ2(1);DQ2(2)];
M3=[DQ3(1);DQ3(2)];
M4=[DQ4(1);DQ4(2)];
M5=[DQ5(1);DQ5(2)];
M6=[DQ6(1);DQ6(2)];
figure;
k=100;
plot(sindqx(1:k),sindqy(1:k),'x');
k=0.8;
axis([-k k -k k],'square');
% sindqx=sindqx*0.8;
% sindqy=sindqy*0.8;
% k=08;hold on,quiver(0,0,sindqx(k),sindqy(k),'k');
% k=26;hold on,quiver(0,0,sindqx(k),sindqy(k),'k');
% k=45;hold on,quiver(0,0,sindqx(k),sindqy(k),'k');
% k=60;hold on,quiver(0,0,sindqx(k),sindqy(k),'k');
% k=76;hold on,quiver(0,0,sindqx(k),sindqy(k),'k');
% k=94;hold on,quiver(0,0,sindqx(k),sindqy(k),'k');
hold on,quiver(0,0,DQ1(1),DQ1(2),'k','linewidth', 2);
hold on,quiver(0,0,DQ2(1),DQ2(2),'k','linewidth', 2);
hold on,quiver(0,0,DQ3(1),DQ3(2),'k','linewidth', 2);
hold on,quiver(0,0,DQ4(1),DQ4(2),'k','linewidth', 2);
hold on,quiver(0,0,DQ5(1),DQ5(2),'k','linewidth', 2);
hold on,quiver(0,0,DQ6(1),DQ6(2),'k','linewidth', 2);
text(DQ1(1),DQ1(2),'V5');
text(DQ2(1),DQ2(2),'V3');
text(DQ3(1),DQ3(2),'V4');
text(DQ4(1),DQ4(2),'V1');
text(DQ5(1),DQ5(2),'V6');
text(DQ6(1),DQ6(2),'V2');
% k=08;text(sindqx(k),sindqy(k),'Vref.1');
% k=26;text(sindqx(k),sindqy(k),'Vref.2');
% k=45;text(sindqx(k),sindqy(k),'Vref.3');
% k=60;text(sindqx(k),sindqy(k),'Vref.4');
% k=76;text(sindqx(k),sindqy(k),'Vref.5');
% k=94;text(sindqx(k),sindqy(k),'Vref.6');
% k=8;hold on,quiver(sindqx(k)-0.04,sindqy(k)-0.04,DQ4(1)-sindqx(k)-0.02,DQ4(2)-sindqy(k)+0.05,'g--');
% k=8;hold on,quiver(sindqx(k)-0.04,sindqy(k)-0.04,DQ6(1)-sindqx(k)+0.02,DQ6(2)-sindqy(k)-0.02,'g--');
% k=8;hold on,quiver(sindqx(k)-0.04,sindqy(k)-0.04,0-sindqx(k)+0.05,0-sindqy(k)+0.035,'g--');
% text(0.46,0.12,'\color{green}Verror.1');
% text(0.37,0.3,'\color{green}Verror.2');
% text(0.2,0.05,'\color{green}Verror.z');


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
Q1=zeros(1,fsw*time);
Q2=zeros(1,fsw*time);
Qz=zeros(1,fsw*time);
D1=zeros(1,fsw*time);
D2=zeros(1,fsw*time);
Dz=zeros(1,fsw*time);
for i=1:1:fsw*time%產生Q1,Q2,Qz
    if mod(sector(i),2)==1%返回1,sector就是奇數；返回0,sector就是偶數
        Q1(i)=((2/3)*cos(angle(i)-(pi/3)*(sector(i)-1))-amp)*T1(i);
        Q2(i)=((2/3)*cos((pi/3)*sector(i)-angle(i))-amp)*T2(i);
        Qz(i)=-amp*Tz(i);
        D1(i)=((2/3)*sin(angle(i)-(pi/3)*(sector(i)-1)))*T1(i);
        D2(i)=((-2/3)*sin((pi/3)*sector(i)-angle(i)))*T2(i);
        Dz(i)=0;
    else
        Q1(i)=((2/3)*cos((pi/3)*sector(i)-angle(i))-amp)*T1(i);
        Q2(i)=((2/3)*cos(angle(i)-(pi/3)*(sector(i)-1))-amp)*T2(i);
        Qz(i)=-amp*Tz(i); 
        D1(i)=((2/3)*sin((pi/3)*sector(i)-angle(i)))*T1(i);
        D2(i)=((-2/3)*sin(angle(i)-(pi/3)*(sector(i)-1)))*T2(i);
        Dz(i)=0;
    end
end
QQQ=Qz+Q1+Q2;
DDD=Dz+D1+D2;
%{
論文圖
i=6;
figure;%0127
y=[0 Qz(i)/2 Qz(i)/2+Q1(i) -Qz(i)/2 0];
x=[0 Tz(i)/2 Tz(i)/2+Tless(i) Tz(i)/2+Tless(i)+Tmore(i) Tz(i)+Tless(i)+Tmore(i)];
subplot(2,1,1),plot(x,y);
y=[0 Dz(i)/2 Dz(i)/2+D1(i) -Dz(i)/2 0];
x=[0 Tz(i)/2 Tz(i)/2+Tless(i) Tz(i)/2+Tless(i)+Tmore(i) Tz(i)+Tless(i)+Tmore(i)];
subplot(2,1,2),plot(x,y);
figure;%0121
y=[0 Qz(i) Qz(i)+Q1(i)/2 -Q1(i)/2 0];
x=[0 Tz(i) Tz(i)+Tless(i)/2 Tz(i)+Tless(i)/2+Tmore(i) Tz(i)+Tless(i)+Tmore(i)];
subplot(2,1,1),plot(x,y);
y=[0 Dz(i) Dz(i)+D1(i)/2 -D1(i)/2 0];
x=[0 Tz(i) Tz(i)+Tless(i)/2 Tz(i)+Tless(i)/2+Tmore(i) Tz(i)+Tless(i)+Tmore(i)];
subplot(2,1,2),plot(x,y);
figure;%7212
y=[0 Qz(i) Qz(i)+Q2(i)/2 -Q2(i)/2 0];
x=[0 Tz(i) Tz(i)+Tmore(i)/2 Tz(i)+Tmore(i)/2+Tless(i) Tz(i)+Tless(i)+Tmore(i)];
subplot(2,1,1),plot(x,y);
y=[0 Dz(i) Dz(i)+D2(i)/2 -D2(i)/2 0];
x=[0 Tz(i) Tz(i)+Tmore(i)/2 Tz(i)+Tmore(i)/2+Tless(i) Tz(i)+Tless(i)+Tmore(i)];
subplot(2,1,2),plot(x,y);
%}
for i=1:fsw*time
    %HDF0127
    P1 = 0;
    P2 = Qz(i)/2;
    P3 = Qz(i)/2+Q1(i);
    P4 = Qz(i)/2+Q1(i)+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i);
    HDFQ1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2;
    HDFQ2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDFQ3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDFQ4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*Tz(i)/2;
    R1 = 0;
    R2 = Dz(i)/2;
    R3 = Dz(i)/2+D1(i);
    R4 = Dz(i)/2+D1(i)+D2(i);
    R5 = Dz(i)+D1(i)+D2(i);
    HDFD1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2;
    HDFD2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i);
    HDFD3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDFD4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*Tz(i)/2;
    %HDF0121
    P1 = 0;
    P2 = Qz(i);
    P3 = Qz(i)+Q1(i)/2;
    P4 = Qz(i)+Q1(i)/2+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i);
    HDF0121Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)/2;
    HDF0121Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDF0121Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T1(i)/2;    
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D1(i)/2;
    R4 = Dz(i)+D1(i)/2+D2(i);
    R5 = Dz(i)+D1(i)+D2(i);
    HDF0121D1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF0121D2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)/2;
    HDF0121D3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDF0121D4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T1(i)/2;
    %HDF7212
    P1 = 0;
    P2 = Qz(i);
    P3 = Qz(i)+Q2(i)/2;
    P4 = Qz(i)+Q2(i)/2+Q1(i);
    P5 = Qz(i)+Q2(i)+Q1(i);
    HDF7212Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF7212Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T2(i)/2;
    HDF7212Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T1(i);
    HDF7212Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)/2;
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D2(i)/2;
    R4 = Dz(i)+D2(i)/2+D1(i);
    R5 = Dz(i)+D2(i)+D1(i);
    HDF7212D1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF7212D2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T2(i)/2;
    HDF7212D3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T1(i);
    HDF7212D4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)/2;
end
HDFQ = HDFQ1+HDFQ2+HDFQ3+HDFQ4;
HDFD = HDFD1+HDFD2+HDFD3+HDFD4;
HDF = HDFQ+HDFD;
HDF0121Q = HDF0121Q1+HDF0121Q2+HDF0121Q3+HDF0121Q4;
HDF0121D = HDF0121D1+HDF0121D2+HDF0121D3+HDF0121D4;
HDF0121 = HDF0121Q+HDF0121D;
HDF7212Q = HDF7212Q1+HDF7212Q2+HDF7212Q3+HDF7212Q4;
HDF7212D = HDF7212D1+HDF7212D2+HDF7212D3+HDF7212D4;
HDF7212 = HDF7212Q+HDF7212D;

threeHDF=[HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
Tzbit=bit/2-T1bit-T2bit;
TTT=T1bit+T2bit+Tzbit;
for i=1:fsw*time%主要的大程式碼
   switch minHDFcase3(i)
       case 1%0127
           out000=repmat([0;0;0],1,floor(Tzbit(i)/2));
           out001=repmat(v1(:,i),1,T1bit(i));
           out011=repmat(v2(:,i),1,T2bit(i));
           out111=repmat([1;1;1],1,Tzbit(i)-floor(Tzbit(i)/2));
           out0127=[out000 out001 out011 out111 out111 out011 out001 out000];
           out(:,bit*(i-1)+1:bit*i)=out0127;
       case 2%0121
           out000=repmat([0;0;0],1,Tzbit(i));
           out001=repmat(v1(:,i),1,floor(T1bit(i)/2));
           out011=repmat(v2(:,i),1,T2bit(i));
           out111=repmat(v1(:,i),1,T1bit(i)-floor(T1bit(i)/2));
           out0121=[out000 out001 out011 out111 out111 out011 out001 out000];
           out(:,bit*(i-1)+1:bit*i)=out0121;
       case 3%7212
           out000=repmat([1;1;1],1,Tzbit(i));
           out001=repmat(v2(:,i),1,floor(T2bit(i)/2));
           out011=repmat(v1(:,i),1,T1bit(i));
           out111=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)/2));
           out7212=[out000 out001 out011 out111 out111 out011 out001 out000];
           out(:,bit*(i-1)+1:bit*i)=out7212;
   end
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

figure,plot(angle*180/pi,Q1);
hold on,plot(angle*180/pi,Q2);
hold on,plot(angle*180/pi,Qz);
figure,plot(angle*180/pi,D1);
hold on,plot(angle*180/pi,D2);
hold on,plot(angle*180/pi,Dz);
figure;%HDF
plot(angle*180/pi,HDF,'r');
hold on,plot(angle*180/pi,HDF0121,'b');
hold on,plot(angle*180/pi,HDF7212,'g');
hold on,plot(angle*180/pi,minHDFvalue3,'m');
title('HDF');
legend('HDF', 'HDF0121','HDF7212')
axis([0 400 0.002 0.02]);

N = time*fsw*bit;
Van=[ 2/3 -1/3 -1/3]*out;
Vbn=[-1/3  2/3 -1/3]*out;
Vcn=[-1/3 -1/3  2/3]*out;
fft_Van = fft(Van,N);
mag_Van = abs(fft_Van)/(length(fft_Van)/2);
fft_Vbn = fft(Vbn,N);
mag_Vbn = abs(fft_Vbn)/(length(fft_Vbn)/2);
fft_Vcn = fft(Vcn,N);
mag_Vcn = abs(fft_Vcn)/(length(fft_Vcn)/2);
w = fsw/N;
k = f*time+1 : f*time : N/2;
Vs1s2s3_Van = mag_Van(k);
Vs1s2s3_Vbn = mag_Vbn(k);
Vs1s2s3_Vcn = mag_Vcn(k);
Vs2s3s4_Van = Vs1s2s3_Van(2:end);
Vs2s3s4_Vbn = Vs1s2s3_Vbn(2:end);
Vs2s3s4_Vcn = Vs1s2s3_Vcn(2:end);
THD_Van = 100*sqrt(sum(Vs2s3s4_Van.^2))/Vs1s2s3_Van(1);
THD_Vbn = 100*sqrt(sum(Vs2s3s4_Vbn.^2))/Vs1s2s3_Vbn(1);
THD_Vcn = 100*sqrt(sum(Vs2s3s4_Vcn.^2))/Vs1s2s3_Vcn(1);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n',THD_Vcn);
Vs2s3s4_Ian = Vs2s3s4_Van./(2:length(Vs2s3s4_Van)+1);
Vs2s3s4_Ibn = Vs2s3s4_Vbn./(2:length(Vs2s3s4_Vbn)+1);
Vs2s3s4_Icn = Vs2s3s4_Vcn./(2:length(Vs2s3s4_Vcn)+1);
ITHD_Van = 100*sqrt(sum(Vs2s3s4_Ian.^2))/Vs1s2s3_Van(1);
ITHD_Vbn = 100*sqrt(sum(Vs2s3s4_Ibn.^2))/Vs1s2s3_Vbn(1);
ITHD_Vcn = 100*sqrt(sum(Vs2s3s4_Icn.^2))/Vs1s2s3_Vcn(1);
fprintf('ITHD1=%f\n',ITHD_Van);
fprintf('ITHD2=%f\n',ITHD_Vbn);
fprintf('ITHD3=%f\n',ITHD_Vcn);
Fdist=sqrt(sum(HDF));
%{
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
% mag_Van_1=[zeros(1,59) mag_Van_1].*36;
% mag_Van_1(60)=18;
% figure,plot(mag_Van_1);
%}
allsin=[s1;s2;s3];