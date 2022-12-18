clear%7/20,7/24ノ程HDF糶
close all
clc
start_time=clock;
time = 0.02;
fsw = 5000;
bit = 10;
amp = 0.5;
f = 50;
n=fsw*time;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t);
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3); 
a=2*pi/3;
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
%plot(sindqy,sindqx);
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

duty=zeros(2,fsw*time);
sector=zeros(1,500);
a0=zeros(1,500);
a7=zeros(1,500);
aless=zeros(1,500);
amore=zeros(1,500);
for i=1:fsw*time%耞sector
if angle(i)>0 && angle(i)<=pi/3;%sector1
    sector(i)=1;
    duty=[M4 M6]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    a0(i)=(1-aless(i)-amore(i))/2;
    a7(i)=a0(i);
end
if angle(i)>pi/3 && angle(i)<=2*pi/3;%sector2
    sector(i)=2;
    duty=[M2 M6]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    a0(i)=(1-aless(i)-amore(i))/2;
    a7(i)=a0(i);
end   
if angle(i)>2*pi/3 && angle(i)<=pi;%sector3
    sector(i)=3;
    duty=[M2 M3]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    a0(i)=(1-aless(i)-amore(i))/2;
    a7(i)=a0(i);
end
if angle(i)>pi && angle(i)<=4*pi/3;%sector4
    sector(i)=4;
    duty=[M1 M3]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    a0(i)=(1-aless(i)-amore(i))/2;
    a7(i)=a0(i);
end
if angle(i)>4*pi/3 && angle(i)<=5*pi/3;%sector5
    sector(i)=5;
    duty=[M1 M5]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    a0(i)=(1-aless(i)-amore(i))/2;
    a7(i)=a0(i);
end 
if angle(i)>5*pi/3 && angle(i)<=2*pi;%sector6
    sector(i)=6;
    duty=[M4 M5]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
    a0(i)=(1-aless(i)-amore(i))/2;
    a7(i)=a0(i);
end
end

Tless=aless;%礚蛤Τ弘非
Tmore=amore;
T0=a0;
T7=a7;
az=a0+a7;
Tz=az;
Q1=zeros(1,fsw*time);
Q2=zeros(1,fsw*time);
Qz=zeros(1,fsw*time);
D1=zeros(1,fsw*time);
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
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(i),2)==1%1,sector碞琌计0,sector碞琌案计
        Q1(i)=((2/3)*cos(angle(i)-(pi/3)*(sector(i)-1))-amp)*Tless(i);
        Q2(i)=((2/3)*cos((pi/3)*sector(i)-angle(i))-amp)*Tmore(i);
        Qz(i)=-amp*Tz(i);
        D1(i)=((2/3)*sin(angle(i)-(pi/3)*(sector(i)-1)))*Tless(i);
        D2(i)=((-2/3)*sin((pi/3)*sector(i)-angle(i)))*Tmore(i);
        Dz(i)=0;
    else
        Q1(i)=((2/3)*cos((pi/3)*sector(i)-angle(i))-amp)*Tless(i);
        Q2(i)=((2/3)*cos(angle(i)-(pi/3)*(sector(i)-1))-amp)*Tmore(i);
        Qz(i)=-amp*Tz(i); 
        D1(i)=((2/3)*sin((pi/3)*sector(i)-angle(i)))*Tless(i);
        D2(i)=((-2/3)*sin(angle(i)-(pi/3)*(sector(i)-1)))*Tmore(i);
        Dz(i)=0;
    end
end
QQQ=Qz+Q1+Q2;
DDD=Dz+D1+D2;
%{
阶ゅ瓜
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
for i=1:1:fsw*time%HDF0127
HDFQ1(i)=((0)^2+(0)*(Qz(i)/2)+(Qz(i)/2)^2)*Tz(i)/2;
HDFQ2(i)=((Qz(i)/2)^2+(Qz(i)/2)*(Qz(i)/2+Q1(i))+(Qz(i)/2+Q1(i))^2)*Tless(i);
HDFQ3(i)=((Qz(i)/2+Q1(i))^2+(Qz(i)/2+Q1(i))*(-Qz(i)/2)+(-Qz(i)/2)^2)*Tmore(i);
HDFQ4(i)=((-Qz(i)/2)^2+(-Qz(i)/2)*(0)+(0)^2)*Tz(i)/2;
HDFQ(i)=HDFQ1(i)+HDFQ2(i)+HDFQ3(i)+HDFQ4(i);
HDFD1(i)=((0)^2+(0)*(Dz(i)/2)+(Dz(i)/2)^2)*Tz(i)/2;
HDFD2(i)=((Dz(i)/2)^2+(Dz(i)/2)*(Dz(i)/2+D1(i))+(Dz(i)/2+D1(i))^2)*Tless(i);
HDFD3(i)=((Dz(i)/2+D1(i))^2+(Dz(i)/2+D1(i))*(-Dz(i)/2)+(-Dz(i)/2)^2)*Tmore(i);
HDFD4(i)=((-Dz(i)/2)^2+(-Dz(i)/2)*(0)+(0)^2)*Tz(i)/2;
HDFD(i)=HDFD1(i)+HDFD2(i)+HDFD3(i)+HDFD4(i);
end
HDF=HDFQ+HDFD;
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
for i=1:1:fsw*time%HDF0121
HDF0121Q1(i)=((0)^2+(0)*(Qz(i))+(Qz(i))^2)*Tz(i);
HDF0121Q2(i)=((Qz(i))^2+(Qz(i))*(Qz(i)+Q1(i)/2)+(Qz(i)+Q1(i)/2)^2)*Tless(i)/2;
HDF0121Q3(i)=((Qz(i)+Q1(i)/2)^2+(Qz(i)+Q1(i)/2)*(-Q1(i)/2)+(-Q1(i)/2)^2)*Tmore(i);
HDF0121Q4(i)=((-Q1(i)/2)^2+(-Q1(i)/2)*(0)+(0)^2)*Tless(i)/2;
HDF0121Q(i)=HDF0121Q1(i)+HDF0121Q2(i)+HDF0121Q3(i)+HDF0121Q4(i);
HDF0121D1(i)=((0)^2+(0)*(Dz(i))+(Dz(i))^2)*Tz(i);
HDF0121D2(i)=((Dz(i))^2+(Dz(i))*(Dz(i)+D1(i)/2)+(Dz(i)+D1(i)/2)^2)*Tless(i)/2;
HDF0121D3(i)=((Dz(i)+D1(i)/2)^2+(Dz(i)+D1(i)/2)*(-D1(i)/2)+(-D1(i)/2)^2)*Tmore(i);
HDF0121D4(i)=((-D1(i)/2)^2+(-D1(i)/2)*(0)+(0)^2)*Tless(i)/2;
HDF0121D(i)=HDF0121D1(i)+HDF0121D2(i)+HDF0121D3(i)+HDF0121D4(i);
end
HDF0121=HDF0121Q+HDF0121D;
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
for i=1:1:fsw*time%HDF7212
HDF7212Q1(i)=((0)^2+(0)*(Qz(i))+(Qz(i))^2)*Tz(i);
HDF7212Q2(i)=((Qz(i))^2+(Qz(i))*(Qz(i)+Q2(i)/2)+(Qz(i)+Q2(i)/2)^2)*Tmore(i)/2;
HDF7212Q3(i)=((Qz(i)+Q2(i)/2)^2+(Qz(i)+Q2(i)/2)*(-Q2(i)/2)+(-Q2(i)/2)^2)*Tless(i);
HDF7212Q4(i)=((-Q2(i)/2)^2+(-Q2(i)/2)*(0)+(0)^2)*Tmore(i)/2;
HDF7212Q(i)=HDF7212Q1(i)+HDF7212Q2(i)+HDF7212Q3(i)+HDF7212Q4(i);
HDF7212D1(i)=((0)^2+(0)*(Dz(i))+(Dz(i))^2)*Tz(i);
HDF7212D2(i)=((Dz(i))^2+(Dz(i))*(Dz(i)+D2(i)/2)+(Dz(i)+D2(i)/2)^2)*Tmore(i)/2;
HDF7212D3(i)=((Dz(i)+D2(i)/2)^2+(Dz(i)+D2(i)/2)*(-D2(i)/2)+(-D2(i)/2)^2)*Tless(i);
HDF7212D4(i)=((-D2(i)/2)^2+(-D2(i)/2)*(0)+(0)^2)*Tmore(i)/2;
HDF7212D(i)=HDF7212D1(i)+HDF7212D2(i)+HDF7212D3(i)+HDF7212D4(i);
end
HDF7212=HDF7212Q+HDF7212D;

partless=round(aless*bit/2);
partmore=round(amore*bit/2);
partz=bit/2-partless-partmore;
ppp=partless+partmore+partz;

threeHDF=[HDF;HDF0121;HDF7212];
minHDFcase=zeros(1,fsw*time);
minHDFvalue=min(threeHDF);
out=zeros(3,fsw*time);
for i=1:1:fsw*time
    for j=1:1:3
        if threeHDF(j,i)==minHDFvalue(i)
            minHDFcase(i)=j;
        end
    end
end
minHDFcase(1,26)=2;
for i=1:fsw*time%璶祘Α絏
   switch minHDFcase(i)
       case 1%0127
           m1=floor(partz/2);
           m2=partless;
           m3=partmore;
           m4=ceil(partz/2);
           m4=2*m4;
    switch sector(i)
        case 1
            for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
            out(:,j)=[0;0;0];
            end
            for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
            out(:,j)=[1;0;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
            out(:,j)=[1;1;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
            out(:,j)=[1;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
            out(:,j)=[1;1;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
            out(:,j)=[1;0;0];
            end
            for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
            out(:,j)=[0;0;0];
            end
         case 2
            for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
            out(:,j)=[0;0;0];
            end
            for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
            out(:,j)=[0;1;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
            out(:,j)=[1;1;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
            out(:,j)=[1;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
            out(:,j)=[1;1;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
            out(:,j)=[0;1;0];
            end
            for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
            out(:,j)=[0;0;0];
            end
         case 3
            for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
            out(:,j)=[0;0;0];
            end
            for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
            out(:,j)=[0;1;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
            out(:,j)=[0;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
            out(:,j)=[1;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
            out(:,j)=[0;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
            out(:,j)=[0;1;0];
            end
            for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
            out(:,j)=[0;0;0];
            end
         case 4
            for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
            out(:,j)=[0;0;0];
            end
            for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
            out(:,j)=[0;0;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
            out(:,j)=[0;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
            out(:,j)=[1;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
            out(:,j)=[0;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
            out(:,j)=[0;0;1];
            end
            for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
            out(:,j)=[0;0;0];
            end
         case 5
            for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
            out(:,j)=[0;0;0];
            end
            for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
            out(:,j)=[0;0;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
            out(:,j)=[1;0;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
            out(:,j)=[1;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
            out(:,j)=[1;0;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
            out(:,j)=[0;0;1];
            end
            for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
            out(:,j)=[0;0;0];
            end
         case 6
            for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
            out(:,j)=[0;0;0];
            end
            for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
            out(:,j)=[1;0;0];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
            out(:,j)=[1;0;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
            out(:,j)=[1;1;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
            out(:,j)=[1;0;1];
            end
            for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
            out(:,j)=[1;0;0];
            end
            for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
            out(:,j)=[0;0;0];
            end
    end
       case 2%0121
           m1=partz;
           m2=floor(partless/2);
           m3=partmore;
           m4=ceil(partless/2);
           m4=2*m4;
       switch sector(i)
           case 1
               for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
               out(:,j)=[0;0;0];
               end
               for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
               out(:,j)=[1;0;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
               out(:,j)=[1;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
               out(:,j)=[1;0;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
               out(:,j)=[1;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
               out(:,j)=[1;0;0];
               end
               for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
               out(:,j)=[0;0;0];
               end
           case 2
               for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
               out(:,j)=[0;0;0];
               end
               for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
               out(:,j)=[0;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
               out(:,j)=[1;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
               out(:,j)=[0;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
               out(:,j)=[1;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
               out(:,j)=[0;1;0];
               end
               for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
               out(:,j)=[0;0;0];
               end
           case 3
               for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
               out(:,j)=[0;0;0];
               end
               for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
               out(:,j)=[0;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
               out(:,j)=[0;1;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
               out(:,j)=[0;1;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
               out(:,j)=[0;1;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
               out(:,j)=[0;1;0];
               end
               for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
               out(:,j)=[0;0;0];
               end
           case 4
               for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
               out(:,j)=[0;0;0];
               end
               for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
               out(:,j)=[0;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
               out(:,j)=[0;1;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
               out(:,j)=[0;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
               out(:,j)=[0;1;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
               out(:,j)=[0;0;1];
               end
               for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
               out(:,j)=[0;0;0];
               end
           case 5
               for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
               out(:,j)=[0;0;0];
               end
               for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
               out(:,j)=[0;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
               out(:,j)=[1;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
               out(:,j)=[0;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
               out(:,j)=[1;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
               out(:,j)=[0;0;1];
               end
               for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
               out(:,j)=[0;0;0];
               end
           case 6
               for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
               out(:,j)=[0;0;0];
               end
               for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
               out(:,j)=[1;0;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
               out(:,j)=[1;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
               out(:,j)=[1;0;0];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
               out(:,j)=[1;0;1];
               end
               for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
               out(:,j)=[1;0;0];
               end
               for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
               out(:,j)=[0;0;0];
               end
       end
       case 3%7212
           m1=partz;
           m2=floor(partmore/2);
           m3=partless;
           m4=ceil(partmore/2);
           m4=2*m4;
           switch sector(i)
               case 1
                   for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
                   out(:,j)=[1;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
                   out(:,j)=[1;1;1];
                   end
               case 2
                   for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
                   out(:,j)=[1;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
                   out(:,j)=[1;1;1];
                   end
               case 3
                   for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
                   out(:,j)=[1;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
                   out(:,j)=[1;1;1];
                   end
               case 4
                   for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
                   out(:,j)=[1;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
                   out(:,j)=[1;1;1];
                   end
               case 5
                   for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
                   out(:,j)=[1;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
                   out(:,j)=[1;1;1];
                   end
               case 6
                   for j=bit*(i-1)+1:1:bit*(i-1)+m1(i);
                   out(:,j)=[1;1;1];
                   end
                   for j=bit*(i-1)+m1(i)+1:1:bit*(i-1)+m1(i)+m2(i);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+1:1:bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+m1(i)+m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+m1(i)+2*m2(i)+2*m3(i)+m4(i)+1:1:bit*(i-1)+bit;
                   out(:,j)=[1;1;1];
                   end
           end

   end
end
% out(2,43)=0;
% out(2,44)=1;
%{
for i=1:fsw*time0121
           switch sector(i)
               case 1
                   for j=bit*(i-1)+1:1:bit*(i-1)+fix(az(i)*bit/2);
                   out(:,j)=[0;0;0];
                   end
                   for j=bit*(i-1)+fix(az(i)*bit/2)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit)+1:1:bit*(i-1)+bit;
                   out(:,j)=[0;0;0];
                   end
               case 2
                   for j=bit*(i-1)+1:1:bit*(i-1)+fix(az(i)*bit/2);
                   out(:,j)=[0;0;0];
                   end
                   for j=bit*(i-1)+fix(az(i)*bit/2)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit);
                   out(:,j)=[1;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit)+1:1:bit*(i-1)+bit;
                   out(:,j)=[0;0;0];
                   end
               case 3
                   for j=bit*(i-1)+1:1:bit*(i-1)+fix(az(i)*bit/2);
                   out(:,j)=[0;0;0];
                   end
                   for j=bit*(i-1)+fix(az(i)*bit/2)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit);
                   out(:,j)=[0;1;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit)+1:1:bit*(i-1)+bit;
                   out(:,j)=[0;0;0];
                   end
               case 4
                   for j=bit*(i-1)+1:1:bit*(i-1)+fix(az(i)*bit/2);
                   out(:,j)=[0;0;0];
                   end
                   for j=bit*(i-1)+fix(az(i)*bit/2)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit);
                   out(:,j)=[0;1;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit)+1:1:bit*(i-1)+bit;
                   out(:,j)=[0;0;0];
                   end
               case 5
                   for j=bit*(i-1)+1:1:bit*(i-1)+fix(az(i)*bit/2);
                   out(:,j)=[0;0;0];
                   end
                   for j=bit*(i-1)+fix(az(i)*bit/2)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit);
                   out(:,j)=[0;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit)+1:1:bit*(i-1)+bit;
                   out(:,j)=[0;0;0];
                   end
               case 6
                   for j=bit*(i-1)+1:1:bit*(i-1)+fix(az(i)*bit/2);
                   out(:,j)=[0;0;0];
                   end
                   for j=bit*(i-1)+fix(az(i)*bit/2)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4)*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i)/2)*bit)+1:1:bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit);
                   out(:,j)=[1;0;1];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+3*aless(i)/4+amore(i))*bit)+1:1:bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit);
                   out(:,j)=[1;0;0];
                   end
                   for j=bit*(i-1)+fix((az(i)/2+aless(i)+amore(i))*bit)+1:1:bit*(i-1)+bit;
                   out(:,j)=[0;0;0];
                   end
           end
end
%}
figure;
out2(1,:)=out(1,2:end);
out2(2,:)=out(2,2:end);
out2(3,:)=out(3,2:end);
out2(1,fsw*time*bit)=0;
out2(2,fsw*time*bit)=0;
out2(3,fsw*time*bit)=0;
out3=[out;out2];
diff1=sum(abs(out(1,:)-out2(1,:)));
diff2=sum(abs(out(2,:)-out2(2,:)));
diff3=sum(abs(out(3,:)-out2(3,:)));
fprintf('out1ち传Ω计=%f\n',diff1);
fprintf('out2ち传Ω计=%f\n',diff2);
fprintf('out3ち传Ω计=%f\n\n',diff3);

subplot(4,1,1),plot(t,s1,t,s2,t,s3);
title('块');
subplot(4,1,2),plot(out(1,:));
title('out1块');
subplot(4,1,3),plot(out(2,:));
title('out2块');
subplot(4,1,4),plot(out(3,:));
title('out3块');

figure;%HDF
plot(angle*180/pi,HDF,'r');
hold on;
plot(angle*180/pi,HDF0121,'b');
hold on;
plot(angle*180/pi,HDF7212,'g');
title('HDF');

N = time*fsw*bit;

fft_Van=fft(out(1,:),N)/N;
mag_Van=abs(fft_Van)*2;

fft_Vbn=fft(out(2,:),N)/N;
mag_Vbn=abs(fft_Vbn)*2;

fft_Vcn=fft(out(3,:),N)/N;
mag_Vcn=abs(fft_Vcn)*2;

i=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
mag_Van_1=mag_Van(1,i);
mag_Vbn_1=mag_Vbn(1,i);
mag_Vcn_1=mag_Vcn(1,i);

THD_Van=100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n\n',THD_Vcn);

n2=2:1:N/2;
mag_Van_2(1,n2-1)=mag_Van(1,n2);
mag_Vbn_2(1,n2-1)=mag_Vbn(1,n2); 
mag_Vcn_2(1,n2-1)=mag_Vcn(1,n2); 

number=size(mag_Van_2);
m2=number(1,2);

for n3=1:1:m2;
mag_Van_3(1,(bit*fsw/N)*n3)=mag_Van_2(1,n3);
mag_Vbn_3(1,(bit*fsw/N)*n3)=mag_Vbn_2(1,n3);
mag_Vcn_3(1,(bit*fsw/N)*n3)=mag_Vcn_2(1,n3);
end

figure
  semilogx(mag_Van_3,'r');
  hold on;
  plot(mag_Van_3,'b');
%   hold on;
%   plot(mag_Vcn_3,'k');


%{
fft_out1 = fft(out(1,:),N); 
mag_out1 = abs(fft_out1)/(length(fft_out1)/2); 
fft_out2 = fft(out(2,:),N); 
mag_out2 = abs(fft_out2)/(length(fft_out2)/2); 
fft_out3 = fft(out(3,:),N); 
mag_out3 = abs(fft_out3)/(length(fft_out3)/2); 
freq = (1:length(fft_out1)/2-1)*fsw*bit/length(fft_out1); 
%freq = (0:150000/2-1)*10; 
draw1=mag_out1(2:1:N/2);%礶fft暗硂妓锣传
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
fprintf('VTHD1=%f\n',THD1);
fprintf('VTHD2=%f\n',THD2);
fprintf('VTHD3=%f\n\n',THD3);

A = [ 2/3 -1/3 -1/3 ; -1/3 2/3 -1/3 ; -1/3 -1/3 2/3];  
B = [ out(1,:) ; out(2,:) ; out(3,:) ];  
C = 48*A*B  ; 
figure;
plot(C(1,1:1:fsw*bit*time));
hold on;
plot(C(2,1:1:fsw*bit*time));
hold on;
plot(C(3,1:1:fsw*bit*time));
title('筿溃');
% figure
% semilogx(C(1,1:1:fsw*bit*time))    
end_time=clock;
execution_time=end_time-start_time;
%}