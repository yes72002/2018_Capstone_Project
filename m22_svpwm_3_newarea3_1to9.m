clear%10/18新3區,新5區--->只有面積
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
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
for i=1:1:n
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
sector=zeros(1,fsw*time);
aless=zeros(1,fsw*time);
amore=zeros(1,fsw*time);
for i=1:fsw*time%判斷sector
if angle(i)>0 && angle(i)<=pi/3;%sector1
    sector(i)=1;
    duty=[M4 M6]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
end
if angle(i)>pi/3 && angle(i)<=2*pi/3;%sector2
    sector(i)=2;
    duty=[M2 M6]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
end   
if angle(i)>2*pi/3 && angle(i)<=pi;%sector3
    sector(i)=3;
    duty=[M2 M3]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
end
if angle(i)>pi && angle(i)<=4*pi/3;%sector4
    sector(i)=4;
    duty=[M1 M3]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
end
if angle(i)>4*pi/3 && angle(i)<=5*pi/3;%sector5
    sector(i)=5;
    duty=[M1 M5]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
end 
if angle(i)>5*pi/3 && angle(i)<=2*pi;%sector6
    sector(i)=6;
    duty=[M4 M5]\[sindqx(i) ; sindqy(i)];
    aless(i)=duty(1);
    amore(i)=duty(2);
end
end
az=1-aless-amore;
a0=az/2;
a7=az/2;
fraction=0.1:0.1:0.9;

Tless=aless;%無限跟有限精準度
Tmore=amore;
Tz=az;
Q1=zeros(1,fsw*time);
Q2=zeros(1,fsw*time);
Qz=zeros(1,fsw*time);
D1=zeros(1,fsw*time);
D2=zeros(1,fsw*time);
Dz=zeros(1,fsw*time);
HDFQ1=zeros(1,fsw*time);
HDFQ2=zeros(1,fsw*time);
HDFQ3=zeros(1,fsw*time);
HDFQ4=zeros(1,fsw*time);
HDFD1=zeros(1,fsw*time);
HDFD2=zeros(1,fsw*time);
HDFD3=zeros(1,fsw*time);
HDFD4=zeros(1,fsw*time);
for i=1:1:fsw*time%產生Q1,Q2,Qz
    if mod(sector(i),2)==1%返回1,sector就是奇數；返回0,sector就是偶數
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
for i=1:1:fsw*time%HDF0127
HDFQ1(i)=((0)^2+(0)*(Qz(i)/2)+(Qz(i)/2)^2)*Tz(i)/2;
HDFQ2(i)=((Qz(i)/2)^2+(Qz(i)/2)*(Qz(i)/2+Q1(i))+(Qz(i)/2+Q1(i))^2)*Tless(i);
HDFQ3(i)=((Qz(i)/2+Q1(i))^2+(Qz(i)/2+Q1(i))*(-Qz(i)/2)+(-Qz(i)/2)^2)*Tmore(i);
HDFQ4(i)=((-Qz(i)/2)^2+(-Qz(i)/2)*(0)+(0)^2)*Tz(i)/2;
HDFD1(i)=((0)^2+(0)*(Dz(i)/2)+(Dz(i)/2)^2)*Tz(i)/2;
HDFD2(i)=((Dz(i)/2)^2+(Dz(i)/2)*(Dz(i)/2+D1(i))+(Dz(i)/2+D1(i))^2)*Tless(i);
HDFD3(i)=((Dz(i)/2+D1(i))^2+(Dz(i)/2+D1(i))*(-Dz(i)/2)+(-Dz(i)/2)^2)*Tmore(i);
HDFD4(i)=((-Dz(i)/2)^2+(-Dz(i)/2)*(0)+(0)^2)*Tz(i)/2;
end
HDFQ=HDFQ1+HDFQ2+HDFQ3+HDFQ4;
HDFD=HDFD1+HDFD2+HDFD3+HDFD4;
HDF=HDFQ+HDFD;
HDF0121Q1=zeros(1,fsw*time);
HDF0121Q2=zeros(9,fsw*time);
HDF0121Q3=zeros(9,fsw*time);
HDF0121Q4=zeros(9,fsw*time);
HDF0121D1=zeros(1,fsw*time);
HDF0121D2=zeros(9,fsw*time);
HDF0121D3=zeros(9,fsw*time);
HDF0121D4=zeros(9,fsw*time);
for k=1:9
for i=1:1:fsw*time%HDF0121
HDF0121Q1(1,i)=((0)^2+(0)*(Qz(i))+(Qz(i))^2)*Tz(i);
HDF0121Q2(k,i)=((Qz(i))^2+(Qz(i))*(Qz(i)+Q1(i)*fraction(k))+(Qz(i)+Q1(i)*fraction(k))^2)*Tless(i)*fraction(k);
HDF0121Q3(k,i)=((Qz(i)+Q1(i)*fraction(k))^2+(Qz(i)+Q1(i)*fraction(k))*(-Q1(i)*(1-fraction(k)))+(-Q1(i)*(1-fraction(k)))^2)*Tmore(i);
HDF0121Q4(k,i)=((-Q1(i)*(1-fraction(k)))^2+(-Q1(i)*(1-fraction(k)))*(0)+(0)^2)*Tless(i)*(1-fraction(k));
HDF0121D1(1,i)=((0)^2+(0)*(Dz(i))+(Dz(i))^2)*Tz(i);
HDF0121D2(k,i)=((Dz(i))^2+(Dz(i))*(Dz(i)+D1(i)*fraction(k))+(Dz(i)+D1(i)*fraction(k))^2)*Tless(i)*fraction(k);
HDF0121D3(k,i)=((Dz(i)+D1(i)*fraction(k))^2+(Dz(i)+D1(i)*fraction(k))*(-D1(i)*(1-fraction(k)))+(-D1(i)*(1-fraction(k)))^2)*Tmore(i);
HDF0121D4(k,i)=((-D1(i)*(1-fraction(k)))^2+(-D1(i)*(1-fraction(k)))*(0)+(0)^2)*Tless(i)*(1-fraction(k));
end
end
for k=1:9
HDF0121Q(k,:)=HDF0121Q1(1,:)+HDF0121Q2(k,:)+HDF0121Q3(k,:)+HDF0121Q4(k,:);
HDF0121D(k,:)=HDF0121D1(1,:)+HDF0121D2(k,:)+HDF0121D3(k,:)+HDF0121D4(k,:);
HDF0121(k,:)=HDF0121Q(k,:)+HDF0121D(k,:);
end
% figure;
% for k=1:9
% plot(angle*180/pi,HDF0121(k,:));
% axis([-inf,inf,0.002,0.02]);
% hold on;
% end

HDF7212Q1=zeros(1,fsw*time);
HDF7212Q2=zeros(9,fsw*time);
HDF7212Q3=zeros(9,fsw*time);
HDF7212Q4=zeros(9,fsw*time);
HDF7212D1=zeros(1,fsw*time);
HDF7212D2=zeros(9,fsw*time);
HDF7212D3=zeros(9,fsw*time);
HDF7212D4=zeros(9,fsw*time);
for k=1:9
for i=1:1:fsw*time%HDF7212
HDF7212Q1(1,i)=((0)^2+(0)*(Qz(i))+(Qz(i))^2)*Tz(i);
HDF7212Q2(k,i)=((Qz(i))^2+(Qz(i))*(Qz(i)+Q2(i)*fraction(k))+(Qz(i)+Q2(i)*fraction(k))^2)*Tmore(i)*fraction(k);
HDF7212Q3(k,i)=((Qz(i)+Q2(i)*fraction(k))^2+(Qz(i)+Q2(i)*fraction(k))*(-Q2(i)*(1-fraction(k)))+(-Q2(i)*(1-fraction(k)))^2)*Tless(i);
HDF7212Q4(k,i)=((-Q2(i)*(1-fraction(k)))^2+(-Q2(i)*(1-fraction(k)))*(0)+(0)^2)*Tmore(i)*(1-fraction(k));
HDF7212D1(1,i)=((0)^2+(0)*(Dz(i))+(Dz(i))^2)*Tz(i);
HDF7212D2(k,i)=((Dz(i))^2+(Dz(i))*(Dz(i)+D2(i)*fraction(k))+(Dz(i)+D2(i)*fraction(k))^2)*Tmore(i)*fraction(k);
HDF7212D3(k,i)=((Dz(i)+D2(i)*fraction(k))^2+(Dz(i)+D2(i)*fraction(k))*(-D2(i)*(1-fraction(k)))+(-D2(i)*(1-fraction(k)))^2)*Tless(i);
HDF7212D4(k,i)=((-D2(i)*(1-fraction(k)))^2+(-D2(i)*(1-fraction(k)))*(0)+(0)^2)*Tmore(i)*(1-fraction(k));
end
end
for k=1:9
HDF7212Q(k,:)=HDF7212Q1(1,:)+HDF7212Q2(k,:)+HDF7212Q3(k,:)+HDF7212Q4(k,:);
HDF7212D(k,:)=HDF7212D1(1,:)+HDF7212D2(k,:)+HDF7212D3(k,:)+HDF7212D4(k,:);
HDF7212(k,:)=HDF7212Q(k,:)+HDF7212D(k,:);
end

threeHDF1=[HDF;HDF0121(1,:);HDF7212(1,:)];
threeHDF2=[HDF;HDF0121(2,:);HDF7212(2,:)];
threeHDF3=[HDF;HDF0121(3,:);HDF7212(3,:)];
threeHDF4=[HDF;HDF0121(4,:);HDF7212(4,:)];
threeHDF5=[HDF;HDF0121(5,:);HDF7212(5,:)];
threeHDF6=[HDF;HDF0121(6,:);HDF7212(6,:)];
threeHDF7=[HDF;HDF0121(7,:);HDF7212(7,:)];
threeHDF8=[HDF;HDF0121(8,:);HDF7212(8,:)];
threeHDF9=[HDF;HDF0121(9,:);HDF7212(9,:)];
minHDFvalue(1,:)=min(threeHDF1);
minHDFvalue(2,:)=min(threeHDF2);
minHDFvalue(3,:)=min(threeHDF3);
minHDFvalue(4,:)=min(threeHDF4);
minHDFvalue(5,:)=min(threeHDF5);
minHDFvalue(6,:)=min(threeHDF6);
minHDFvalue(7,:)=min(threeHDF7);
minHDFvalue(8,:)=min(threeHDF8);
minHDFvalue(9,:)=min(threeHDF9);
for k=1:9
area(k,:)=HDF(1,:)-minHDFvalue(k,:);
allarea(k)=sum(area(k,:));
end
fprintf('allarea=%f\n',allarea);
point=zeros(1,9);
for k=1:9
for i=1:n%1
    if k==1
        threeHDF=threeHDF1;
    elseif k==2
            threeHDF=threeHDF2;
    elseif k==3
            threeHDF=threeHDF3;
    elseif k==4
            threeHDF=threeHDF4;
    elseif k==5
            threeHDF=threeHDF5;
    elseif k==6
            threeHDF=threeHDF6;
    elseif k==7
            threeHDF=threeHDF7;
    elseif k==8
            threeHDF=threeHDF8;
    elseif k==9
            threeHDF=threeHDF9;
    end
    if threeHDF(1,i)==minHDFvalue(k,i)
        minHDFcase(k,i)=1;
    elseif threeHDF(2,i)==minHDFvalue(k,i)
        minHDFcase(k,i)=2;
        point(k)=point(k)+1;
    elseif threeHDF(3,i)==minHDFvalue(k,i)
        minHDFcase(k,i)=3;
        point(k)=point(k)+1;
    end
end
end
fprintf('point=%f\n',point);

% figure;
% plot(angle*180/pi,minHDFvalue,'r');
% axis([-inf,inf,0.002,0.02]);
% hold on;
% plot(angle*180/pi,HDF,'b');
% figure;
% subplot(4,1,1),plot(t,s1,t,s2,t,s3);
% title('三相輸入');
% figure;%HDF
% plot(angle*180/pi,HDF,'r');
% hold on;
% plot(angle*180/pi,HDF0121,'b');
% hold on;
% plot(angle*180/pi,HDF7212,'g');
% title('HDF');
% axis([-inf,inf,0.002,0.02]);