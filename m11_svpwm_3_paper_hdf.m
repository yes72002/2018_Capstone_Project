clear%7/19誤差HDF試試論文
close all
clc
time = 0.1;
fsw = 5000;
bit = 300;
amp = 0.5;
f = 60;
t = 0:1/fsw:time-1/fsw; 
T=1/fsw;
s1 = amp*sin(2*pi*f*t); 
s2 = amp*sin(2*pi*f*t+2*pi/3); 
s3 = amp*sin(2*pi*f*t+4*pi/3); 
subplot(4,1,1),plot(t,s1,t,s2,t,s3);

mapping=2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1 1 1];
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
out=zeros(3,fsw*time);
out1=zeros(1,fsw*bit*time);
out2=zeros(1,fsw*bit*time);
out3=zeros(1,fsw*bit*time);
for i=1:fsw*time%判斷sector並生成a0,aless,amore,a7,sector
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
if angle(i)>pi && angle(i)<=4*pi/3;%sector4
    sector(i)=4;
    duty=[M1 M3]\[sindqx(i) ; sindqy(i)];
    a1=duty(1);
    a3=duty(2);
    a0(i)=(1-a1-a3)/2;
    a7(i)=a0(i);
    aless(i)=a1;
    amore(i)=a3;
end
if angle(i)>4*pi/3 && angle(i)<=5*pi/3;%sector5
    sector(i)=5;
    duty=[M1 M5]\[sindqx(i) ; sindqy(i)];
    a1=duty(1);
    a5=duty(2);
    a0(i)=(1-a1-a5)/2;
    a7(i)=a0(i);
    aless(i)=a1;
    amore(i)=a5;
end 
if angle(i)>5*pi/3 && angle(i)<=2*pi;%sector6
    sector(i)=6;
    duty=[M4 M5]\[sindqx(i) ; sindqy(i)];
    a4=duty(1);
    a5=duty(2);
    a0(i)=(1-a4-a5)/2;
    a7(i)=a0(i);
    aless(i)=a4;
    amore(i)=a5;
end
end
for i=1:fsw*time%產生方波
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

figure;
% plot(sindqx,sindqy);
% plot(a0);
% plot(angle);
M=[M1 M2 M3 M4 M5 M6];
Mx=M(1,1:1:6);
My=M(2,1:1:6);
quiver(zeros(1,6),zeros(1,6),Mx,My);
hold on;
v=20;
quiver(zeros(1,v),zeros(1,v),sindqx(1:1:v),sindqy(1:1:v));

Tless=T*aless;%無限&有限精準度
Tmore=T*amore;
T0=T*a0;
T7=T*a7;
VerrorD=zeros(1,fsw*time);
VerrorD2=zeros(1,fsw*time);
VerrorDz=zeros(1,fsw*time);
VerrorQ1=zeros(1,fsw*time);
VerrorQ2=zeros(1,fsw*time);
VerrorQz=zeros(1,fsw*time);
Q1=zeros(1,fsw*time);
Q2=zeros(1,fsw*time);
Qz=zeros(1,fsw*time);
D=zeros(1,fsw*time);
D2=zeros(1,fsw*time);
Dz=zeros(1,fsw*time);
HDFQ=zeros(1,fsw*time);
HDFQ1=zeros(1,fsw*time);
HDFQ2=zeros(1,fsw*time);
HDFQ3=zeros(1,fsw*time);
HDFQ4=zeros(1,fsw*time);
HDFD=zeros(1,fsw*time);
HDFD1=zeros(1,fsw*time);
HDFD2=zeros(1,fsw*time);
HDFD3=zeros(1,fsw*time);
HDFD4=zeros(1,fsw*time);
for i=1:1:fsw*time
    if mod(sector(i),2)==1%返回1,sector就是奇數；返回0,sector就是偶數
        VerrorD(i)=(2/3)*sin(angle(i)-(pi/3)*(sector(i)-1));
        VerrorD2(i)=(-2/3)*sin((pi/3)*sector(i)-angle(i));
        VerrorDz(i)=0;
        VerrorQ1(i)=(2/3)*cos(angle(i)-(pi/3)*(sector(i)-1))-amp;
        VerrorQ2(i)=(2/3)*cos((pi/3)*sector(i)-angle(i))-amp;
        VerrorQz(i)=-amp; 
    else
        VerrorD(i)=(2/3)*sin((pi/3)*sector(i)-angle(i));
        VerrorD2(i)=(-2/3)*sin(angle(i)-(pi/3)*(sector(i)-1));
        VerrorDz(i)=0;
        VerrorQ1(i)=(2/3)*cos((pi/3)*sector(i)-angle(i))-amp;
        VerrorQ2(i)=(2/3)*cos(angle(i)-(pi/3)*(sector(i)-1))-amp;
        VerrorQz(i)=-amp; 
    end
    Q1(i)=VerrorQ1(i)*Tless(i);
    Q2(i)=VerrorQ2(i)*Tmore(i);
    Qz(i)=VerrorQz(i)*(T0(i)+T7(i));
    D(i)=VerrorD(i)*Tless(i);
end;
% figure;
% plot(VerrorQ1);
% figure;
% plot(VerrorQ2);
% figure;
% QQQ=Q1+Q2+Qz;
% plot(QQQ);
% DDD=D1+D2+Dz;
% figure;
% plot(DDD);
% figure;
% plot(VerrorQ1);
% hold on;
% % figure;
% plot(VerrorQ2);
% hold on;
% % figure;
% plot(VerrorQz);
% hold on;
% figure;
% plot(VerrorD1);
% hold on;
% figure;
% plot(VerrorD2);
% hold on;
% figure;
% plot(VerrorDz);
for i=1:1:fsw*time
HDFQ1(i)=((0)^2+(0)*(Qz(i)/2)+(Qz(i)/2)^2)*(T0(i)+T7(i))/2;
HDFQ2(i)=((Qz(i)/2)^2+(Qz(i)/2)*(Qz(i)/2+Q1(i))+(Qz(i)/2+Q1(i))^2)*Tless(i);
HDFQ3(i)=((Qz(i)/2+Q1(i))^2+(Qz(i)/2+Q1(i))*(-Qz(i)/2)+(-Qz(i)/2)^2)*Tmore(i);
HDFQ4(i)=((-Qz(i)/2)^2+(-Qz(i)/2)*(0)+(0)^2)*(T0(i)+T7(i))/2;
HDFQ(i)=HDFQ1(i)+HDFQ2(i)+HDFQ3(i)+HDFQ4(i);
end
for i=1:1:fsw*time
HDFD1(i)=((0)^2+(0)*(0)+(0)^2)*(T0(i)+T7(i))/2;
HDFD2(i)=((0)^2+(0)*(D(i))+(D(i))^2)*Tless(i);
HDFD3(i)=((D(i))^2+(D(i))*(0)+(0)^2)*Tmore(i);
HDFD4(i)=((0)^2+(0)*(0)+(0)^2)*(T0(i)+T7(i))/2;
HDFD(i)=HDFD1(i)+HDFD2(i)+HDFD3(i)+HDFD4(i);
end
HDF=HDFQ+HDFD;
figure;
plot(angle,HDF,'r');
hold on;

HDF0121Q=zeros(1,fsw*time);
HDF0121Q1=zeros(1,fsw*time);
HDF0121Q2=zeros(1,fsw*time);
HDF0121Q3=zeros(1,fsw*time);
HDF0121Q4=zeros(1,fsw*time);
HDF0121D=zeros(1,fsw*time);
HDF0121D1=zeros(1,fsw*time);
HDF0121D2=zeros(1,fsw*time);
HDF0121D3=zeros(1,fsw*time);
HDF0121D4=zeros(1,fsw*time);
for i=1:1:fsw*time
HDF0121Q1(i)=((0)^2+(0)*(Qz(i))+(Qz(i))^2)*(T0(i)+T7(i));
HDF0121Q2(i)=((Qz(i))^2+(Qz(i))*(Qz(i)+Q1(i)/2)+(Qz(i)+Q1(i)/2)^2)*Tless(i)/2;
HDF0121Q3(i)=((Qz(i)+Q1(i)/2)^2+(Qz(i)+Q1(i)/2)*(-Q1(i)/2)+(-Q1(i)/2)^2)*Tmore(i);
HDF0121Q4(i)=((-Q1(i)/2)^2+(-Q1(i)/2)*(0)+(0)^2)*Tless(i)/2;
HDF0121Q(i)=HDF0121Q1(i)+HDF0121Q2(i)+HDF0121Q3(i)+HDF0121Q4(i);
end
for i=1:1:fsw*time
HDF0121D1(i)=((0)^2+(0)*(0)+(0)^2)*(T0(i)+T7(i));
HDF0121D2(i)=((0)^2+(0)*(D(i)/2)+(D(i)/2)^2)*Tless(i)/2;
HDF0121D3(i)=((D(i)/2)^2+(D(i)/2)*(-D(i)/2)+(-D(i)/2)^2)*Tmore(i);
HDF0121D4(i)=((-D(i)/2)^2+(-D(i)/2)*(0)+(0)^2)*Tless(i)/2;
HDF0121D(i)=HDF0121D1(i)+HDF0121D2(i)+HDF0121D3(i)+HDF0121D4(i);
end
HDF0121=HDF0121Q+HDF0121D;
plot(angle,HDF0121,'b');
hold on;

HDF7212Q=zeros(1,fsw*time);
HDF7212Q1=zeros(1,fsw*time);
HDF7212Q2=zeros(1,fsw*time);
HDF7212Q3=zeros(1,fsw*time);
HDF7212Q4=zeros(1,fsw*time);
HDF7212D=zeros(1,fsw*time);
HDF7212D1=zeros(1,fsw*time);
HDF7212D2=zeros(1,fsw*time);
HDF7212D3=zeros(1,fsw*time);
HDF7212D4=zeros(1,fsw*time);
for i=1:1:fsw*time
HDF7212Q1(i)=((0)^2+(0)*(Qz(i))+(Qz(i))^2)*(T0(i)+T7(i));
HDF7212Q2(i)=((Qz(i))^2+(Qz(i))*(Qz(i)+Q2(i)/2)+(Qz(i)+Q2(i)/2)^2)*Tmore(i)/2;
HDF7212Q3(i)=((Qz(i)+Q2(i)/2)^2+(Qz(i)+Q2(i)/2)*(-Q2(i)/2)+(-Q2(i)/2)^2)*Tless(i);
HDF7212Q4(i)=((-Q2(i)/2)^2+(-Q2(i)/2)*(0)+(0)^2)*Tmore(i)/2;
HDF7212Q(i)=HDF7212Q1(i)+HDF7212Q2(i)+HDF7212Q3(i)+HDF7212Q4(i);
end
for i=1:1:fsw*time
HDF7212D1(i)=((0)^2+(0)*(0)+(0)^2)*(T0(i)+T7(i));
HDF7212D2(i)=((0)^2+(0)*(-D(i)/2)+(-D(i)/2)^2)*Tmore(i)/2;
HDF7212D3(i)=((-D(i)/2)^2+(-D(i)/2)*(D(i)/2)+(D(i)/2)^2)*Tless(i);
HDF7212D4(i)=((D(i)/2)^2+(D(i)/2)*(0)+(0)^2)*Tmore(i)/2;
HDF7212D(i)=HDF7212D1(i)+HDF7212D2(i)+HDF7212D3(i)+HDF7212D4(i);
end
HDF7212=HDF7212Q+HDF7212D;
plot(angle,HDF7212,'g');

N = time*fsw*bit;
fft_out1 = fft(out1,N); 
mag_out1 = abs(fft_out1)/(length(fft_out1)/2); 
fft_out2 = fft(out2,N); 
mag_out2 = abs(fft_out2)/(length(fft_out2)/2); 
fft_out3 = fft(out3,N); 
mag_out3 = abs(fft_out3)/(length(fft_out3)/2); 
freq = (0:length(fft_out1)/2-1)*fsw*bit/length(fft_out1); 
%freq = (0:150000/2-1)*10; 
draw1=mag_out1(1:1:N/2);%為了畫fft才做這樣的轉換
draw2=mag_out2(1:1:N/2);
draw3=mag_out3(1:1:N/2);
for i=2:N/2
  freq(i-1)=freq(i);
  draw1(i-1)=draw1(i);
  draw2(i-1)=draw2(i);
  draw2(i-1)=draw2(i);
end
figure;
% subplot(3,1,1),plot(freq,draw1);
% subplot(3,1,2),plot(freq,draw2);
% subplot(3,1,3),plot(freq,draw3);
plot(freq,draw1);
figure;
semilogx(mag_out1(1:end/2));
k=f/10;
harm1 = mag_out1( 2*k+1: k: end/2 ); 
THD1 = 100*sqrt(  sum( harm1.^2)/mag_out1(k+1)^2 );
harm2 = mag_out2( 2*k+1: k: end/2 ); 
THD2 = 100*sqrt(  sum( harm2.^2)/mag_out2(k+1)^2 );
harm3 = mag_out3( 2*k+1: k: end/2 ); 
THD3 = 100*sqrt(  sum( harm3.^2)/mag_out3(k+1)^2 );

A = [ 2/3 -1/3 -1/3 ; -1/3 2/3 -1/3 ; -1/3 -1/3 2/3];%相電壓  
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
