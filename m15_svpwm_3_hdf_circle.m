clear%7/20,7/24,7/26利用最低的HDF寫出序列,圓的分布
close all
clc
time = 0.02;
fsw = 15000;
bit = 400;
amp = 0.59;%0.59
f = 50;
t = 0:1/fsw:time-1/fsw; 
interval=0.01;
r=fix(amp/interval);
s1=zeros(r,fsw*time);
s2=zeros(r,fsw*time);
s3=zeros(r,fsw*time);
sindq=zeros(r,fsw*time);
sindqx=zeros(r,fsw*time);
sindqy=zeros(r,fsw*time);
for j=1:1:r
   amp=interval*(j-1);
   s1(j,:) = amp*sin(2*pi*f*t);
   s2(j,:) = amp*sin(2*pi*f*t+2*pi/3);
   s3(j,:) = amp*sin(2*pi*f*t+4*pi/3); 
end
mapping=2/3*[1 -1/2 -1/2;0 sqrt(3)/2 -sqrt(3)/2;1 1 1];
for j=1:1:r
sindq(3*j-2,:)=mapping(1,:)*[s1(j,:);s2(j,:);s3(j,:)];
sindq(3*j-1,:)=mapping(2,:)*[s1(j,:);s2(j,:);s3(j,:)];
sindq(3*j,:)=mapping(3,:)*[s1(j,:);s2(j,:);s3(j,:)];
sindqx(j,:)=sindq(3*j-2,:);
sindqy(j,:)=sindq(3*j-1,:);
end
angle=atan2(sindqy,sindqx);
for i=1:1:fsw*time
    if angle(2,i)<0
        angle(2,i)=angle(2,i)+2*pi;
    end
end
for j=3:1:r
angle(j,:)=angle(2,:);
end
M4=[2/3 ; 0];
M6=[1/3 ; sqrt(3)/3];
M2=[-1/3 ; sqrt(3)/3];
M3=[-2/3 ; 0];
M1=[-1/3 ; -sqrt(3)/3];
M5=[1/3 ; -sqrt(3)/3];

duty=zeros(2*r,fsw*time);
sector=zeros(r,fsw*time);
a0=zeros(r,fsw*time);
a7=zeros(r,fsw*time);
aless=zeros(r,fsw*time);
amore=zeros(r,fsw*time);
for j=1:1:r
for i=1:fsw*time%判斷sector
if angle(j,i)>0 && angle(j,i)<=pi/3;%sector1
    sector(j,i)=1;
    duty=[M4 M6]\[sindqx(j,i) ; sindqy(j,i)];
    a4=duty(1);
    a6=duty(2);
    a0(j,i)=(1-a4-a6)/2;
    a7(j,i)=a0(j,i);
    aless(j,i)=a4;
    amore(j,i)=a6;
end
if angle(j,i)>pi/3 && angle(j,i)<=2*pi/3;%sector2
    sector(j,i)=2;
    duty=[M2 M6]\[sindqx(j,i) ; sindqy(j,i)];
    a2=duty(1);
    a6=duty(2);
    a0(j,i)=(1-a2-a6)/2;
    a7(j,i)=a0(j,i);
    aless(j,i)=a2;
    amore(j,i)=a6;
end   
if angle(j,i)>2*pi/3 && angle(j,i)<=pi;%sector3
    sector(j,i)=3;
    duty=[M2 M3]\[sindqx(j,i) ; sindqy(j,i)];
    a2=duty(1);
    a3=duty(2);
    a0(j,i)=(1-a2-a3)/2;
    a7(j,i)=a0(j,i);
    aless(j,i)=a2;
    amore(j,i)=a3;
end
if angle(j,i)>pi && angle(j,i)<=4*pi/3;%sector4
    sector(j,i)=4;
    duty=[M1 M3]\[sindqx(j,i) ; sindqy(j,i)];
    a1=duty(1);
    a3=duty(2);
    a0(j,i)=(1-a1-a3)/2;
    a7(j,i)=a0(j,i);
    aless(j,i)=a1;
    amore(j,i)=a3;
end
if angle(j,i)>4*pi/3 && angle(j,i)<=5*pi/3;%sector5
    sector(j,i)=5;
    duty=[M1 M5]\[sindqx(j,i) ; sindqy(j,i)];
    a1=duty(1);
    a5=duty(2);
    a0(j,i)=(1-a1-a5)/2;
    a7(j,i)=a0(j,i);
    aless(j,i)=a1;
    amore(j,i)=a5;
end 
if angle(j,i)>5*pi/3 && angle(j,i)<=2*pi;%sector6
    sector(j,i)=6;
    duty=[M4 M5]\[sindqx(j,i) ; sindqy(j,i)];
    a4=duty(1);
    a5=duty(2);
    a0(j,i)=(1-a4-a5)/2;
    a7(j,i)=a0(j,i);
    aless(j,i)=a4;
    amore(j,i)=a5;
end
end
end

Tless=aless;%無限跟有限精準度
Tmore=amore;
T0=a0;
T7=a7;
az=a0+a7;
Tz=az;
Q1=zeros(r,fsw*time);
Q2=zeros(r,fsw*time);
Qz=zeros(r,fsw*time);
D1=zeros(r,fsw*time);
D2=zeros(r,fsw*time);
Dz=zeros(r,fsw*time);
HDFQ1=zeros(r,fsw*time);
HDFQ2=zeros(r,fsw*time);
HDFQ3=zeros(r,fsw*time);
HDFQ4=zeros(r,fsw*time);
HDFD1=zeros(r,fsw*time);
HDFD2=zeros(r,fsw*time);
HDFD3=zeros(r,fsw*time);
HDFD4=zeros(r,fsw*time);
for j=1:1:r
    amp=interval*(j-1);
for i=1:1:fsw*time%產生Q1,Q2,Qz
    if mod(sector(j,i),2)==1%返回1,sector就是奇數；返回0,sector就是偶數
        Q1(j,i)=((2/3)*cos(angle(j,i)-(pi/3)*(sector(j,i)-1))-amp)*Tless(j,i);
        Q2(j,i)=((2/3)*cos((pi/3)*sector(j,i)-angle(j,i))-amp)*Tmore(j,i);
        Qz(j,i)=-amp*Tz(j,i);
        D1(j,i)=((2/3)*sin(angle(j,i)-(pi/3)*(sector(j,i)-1)))*Tless(j,i);
        D2(j,i)=((-2/3)*sin((pi/3)*sector(j,i)-angle(j,i)))*Tmore(j,i);
        Dz(j,i)=0;
    else
        Q1(j,i)=((2/3)*cos((pi/3)*sector(j,i)-angle(j,i))-amp)*Tless(j,i);
        Q2(j,i)=((2/3)*cos(angle(j,i)-(pi/3)*(sector(j,i)-1))-amp)*Tmore(j,i);
        Qz(j,i)=-amp*Tz(j,i); 
        D1(j,i)=((2/3)*sin((pi/3)*sector(j,i)-angle(j,i)))*Tless(j,i);
        D2(j,i)=((-2/3)*sin(angle(j,i)-(pi/3)*(sector(j,i)-1)))*Tmore(j,i);
        Dz(j,i)=0;
    end
end
end
QQQ=Q1+Q2+Qz;
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
for j=1:1:r
for i=1:1:fsw*time%HDF0127
HDFQ1(j,i)=((0)^2+(0)*(Qz(j,i)/2)+(Qz(j,i)/2)^2)*Tz(j,i)/2;
HDFQ2(j,i)=((Qz(j,i)/2)^2+(Qz(j,i)/2)*(Qz(j,i)/2+Q1(j,i))+(Qz(j,i)/2+Q1(j,i))^2)*Tless(j,i);
HDFQ3(j,i)=((Qz(j,i)/2+Q1(j,i))^2+(Qz(j,i)/2+Q1(j,i))*(-Qz(j,i)/2)+(-Qz(j,i)/2)^2)*Tmore(j,i);
HDFQ4(j,i)=((-Qz(j,i)/2)^2+(-Qz(j,i)/2)*(0)+(0)^2)*Tz(j,i)/2;
HDFD1(j,i)=((0)^2+(0)*(Dz(j,i)/2)+(Dz(j,i)/2)^2)*Tz(j,i)/2;
HDFD2(j,i)=((Dz(j,i)/2)^2+(Dz(j,i)/2)*(Dz(j,i)/2+D1(j,i))+(Dz(j,i)/2+D1(j,i))^2)*Tless(j,i);
HDFD3(j,i)=((Dz(j,i)/2+D1(j,i))^2+(Dz(j,i)/2+D1(j,i))*(-Dz(j,i)/2)+(-Dz(j,i)/2)^2)*Tmore(j,i);
HDFD4(j,i)=((-Dz(j,i)/2)^2+(-Dz(j,i)/2)*(0)+(0)^2)*Tz(j,i)/2;
end
end
HDFQ=HDFQ1+HDFQ2+HDFQ3+HDFQ4;
HDFD=HDFD1+HDFD2+HDFD3+HDFD4;
HDF=HDFQ+HDFD;
HDF0121Q1=zeros(r,fsw*time);
HDF0121Q2=zeros(r,fsw*time);
HDF0121Q3=zeros(r,fsw*time);
HDF0121Q4=zeros(r,fsw*time);
HDF0121D1=zeros(r,fsw*time);
HDF0121D2=zeros(r,fsw*time);
HDF0121D3=zeros(r,fsw*time);
HDF0121D4=zeros(r,fsw*time);
for j=1:1:r
for i=1:1:fsw*time%HDF0121
HDF0121Q1(j,i)=((0)^2+(0)*(Qz(j,i))+(Qz(j,i))^2)*Tz(j,i);
HDF0121Q2(j,i)=((Qz(j,i))^2+(Qz(j,i))*(Qz(j,i)+Q1(j,i)/2)+(Qz(j,i)+Q1(j,i)/2)^2)*Tless(j,i)/2;
HDF0121Q3(j,i)=((Qz(j,i)+Q1(j,i)/2)^2+(Qz(j,i)+Q1(j,i)/2)*(-Q1(j,i)/2)+(-Q1(j,i)/2)^2)*Tmore(j,i);
HDF0121Q4(j,i)=((-Q1(j,i)/2)^2+(-Q1(j,i)/2)*(0)+(0)^2)*Tless(j,i)/2;
HDF0121D1(j,i)=((0)^2+(0)*(Dz(j,i))+(Dz(j,i))^2)*Tz(j,i);
HDF0121D2(j,i)=((Dz(j,i))^2+(Dz(j,i))*(Dz(j,i)+D1(j,i)/2)+(Dz(j,i)+D1(j,i)/2)^2)*Tless(j,i)/2;
HDF0121D3(j,i)=((Dz(j,i)+D1(j,i)/2)^2+(Dz(j,i)+D1(j,i)/2)*(-D1(j,i)/2)+(-D1(j,i)/2)^2)*Tmore(j,i);
HDF0121D4(j,i)=((-D1(j,i)/2)^2+(-D1(j,i)/2)*(0)+(0)^2)*Tless(j,i)/2;
end
end
HDF0121Q=HDF0121Q1+HDF0121Q2+HDF0121Q3+HDF0121Q4;
HDF0121D=HDF0121D1+HDF0121D2+HDF0121D3+HDF0121D4;
HDF0121=HDF0121Q+HDF0121D;
HDF7212Q1=zeros(r,fsw*time);
HDF7212Q2=zeros(r,fsw*time);
HDF7212Q3=zeros(r,fsw*time);
HDF7212Q4=zeros(r,fsw*time);
HDF7212D1=zeros(r,fsw*time);
HDF7212D2=zeros(r,fsw*time);
HDF7212D3=zeros(r,fsw*time);
HDF7212D4=zeros(r,fsw*time);
for j=1:1:r
for i=1:1:fsw*time%HDF7212
HDF7212Q1(j,i)=((0)^2+(0)*(Qz(j,i))+(Qz(j,i))^2)*Tz(j,i);
HDF7212Q2(j,i)=((Qz(j,i))^2+(Qz(j,i))*(Qz(j,i)+Q2(j,i)/2)+(Qz(j,i)+Q2(j,i)/2)^2)*Tmore(j,i)/2;
HDF7212Q3(j,i)=((Qz(j,i)+Q2(j,i)/2)^2+(Qz(j,i)+Q2(j,i)/2)*(-Q2(j,i)/2)+(-Q2(j,i)/2)^2)*Tless(j,i);
HDF7212Q4(j,i)=((-Q2(j,i)/2)^2+(-Q2(j,i)/2)*(0)+(0)^2)*Tmore(j,i)/2;
HDF7212D1(j,i)=((0)^2+(0)*(Dz(j,i))+(Dz(j,i))^2)*Tz(j,i);
HDF7212D2(j,i)=((Dz(j,i))^2+(Dz(j,i))*(Dz(j,i)+D2(j,i)/2)+(Dz(j,i)+D2(j,i)/2)^2)*Tmore(j,i)/2;
HDF7212D3(j,i)=((Dz(j,i)+D2(j,i)/2)^2+(Dz(j,i)+D2(j,i)/2)*(-D2(j,i)/2)+(-D2(j,i)/2)^2)*Tless(j,i);
HDF7212D4(j,i)=((-D2(j,i)/2)^2+(-D2(j,i)/2)*(0)+(0)^2)*Tmore(j,i)/2;
end
end
HDF7212Q=HDF7212Q1+HDF7212Q2+HDF7212Q3+HDF7212Q4;
HDF7212D=HDF7212D1+HDF7212D2+HDF7212D3+HDF7212D4;
HDF7212=HDF7212Q+HDF7212D;
threeHDF=zeros(3*r,fsw*time);
minHDFcase=zeros(r,fsw*time);
minHDFvalue=zeros(r,fsw*time);
for j=1:1:r
    threeHDF(3*j-2,:)=HDF(j,:);
    threeHDF(3*j-1,:)=HDF0121(j,:);
    threeHDF(3*j,:)=HDF7212(j,:);
    minHDFvalue(j,:)=min(threeHDF(3*j-2:3*j,:));
end
%out=zeros(3,fsw*time);
for j=1:1:r
    for i=1:1:fsw*time
        if threeHDF(3*j-2,i)==minHDFvalue(j,i)
            minHDFcase(j,i)=1;
        elseif threeHDF(3*j-1,i)==minHDFvalue(j,i)
            minHDFcase(j,i)=2;
        elseif threeHDF(3*j,i)==minHDFvalue(j,i)
            minHDFcase(j,i)=3;
        end
    end
end


% figure;%HDF
% plot(angle*180/pi,HDF,'r');
% hold on;
% plot(angle*180/pi,HDF0121,'b');
% hold on;
% plot(angle*180/pi,HDF7212,'g');
sinHDFdistribute1x=zeros(r,fsw*time);
sinHDFdistribute1y=zeros(r,fsw*time);
sinHDFdistribute2x=zeros(r,fsw*time);
sinHDFdistribute2y=zeros(r,fsw*time);
sinHDFdistribute3x=zeros(r,fsw*time);
sinHDFdistribute3y=zeros(r,fsw*time);
figure;
for j=1:1:r
    plot(sindqx(j,:),sindqy(j,:),'b.');
    hold on;
end
 title('sinwave');
 axis([-0.6 0.6 -0.6 0.6],'square')
figure;%HDF圓的分布
for j=1:1:r
for i=1:1:fsw*time
    switch minHDFcase(j,i)
        case 1
        sinHDFdistribute1x(j,i)=sindqx(j,i);
        sinHDFdistribute1y(j,i)=sindqy(j,i);
        case 2
        sinHDFdistribute2x(j,i)=sindqx(j,i);
        sinHDFdistribute2y(j,i)=sindqy(j,i);
        case 3
        sinHDFdistribute3x(j,i)=sindqx(j,i);
        sinHDFdistribute3y(j,i)=sindqy(j,i);
    end
end
end
for j=1:1:r
plot(sinHDFdistribute1x(j,:),sinHDFdistribute1y(j,:),'r.');
hold on;
plot(sinHDFdistribute2x(j,:),sinHDFdistribute2y(j,:),'b.');
hold on;
plot(sinHDFdistribute3x(j,:),sinHDFdistribute3y(j,:),'k.');
hold on;
axis([-0.6 0.6 -0.6 0.6],'square')
legend('HDF', 'HDF0121','HDF7212')
end
title('hdf circle');
 
Fdist=sqrt(sum(HDF'));
Fdist0121=sqrt(sum(HDF0121'));
Fdist7212=sqrt(sum(HDF7212'));
figure,plot(Fdist,'r');
hold on,plot(Fdist0121,'b');
hold on,plot(Fdist7212,'g');
title('Fdist');
