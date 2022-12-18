clear%11/19穝跋
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.5;
fraction3=0.5;
f = 50;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t);
s2 = amp*cos(2*pi*f*t-2*pi/7); 
s3 = amp*cos(2*pi*f*t-4*pi/7);
s4 = amp*cos(2*pi*f*t-6*pi/7); 
s5 = amp*cos(2*pi*f*t-8*pi/7); 
s6 = amp*cos(2*pi*f*t-10*pi/7); 
s7 = amp*cos(2*pi*f*t-12*pi/7); 
a = 2*pi/7;
mapping=2/7*[...
    cos(0) cos(a) cos(2*a) cos(3*a) cos(4*a) cos(5*a) cos(6*a);...
    sin(0) sin(a) sin(2*a) sin(3*a) sin(4*a) sin(5*a) sin(6*a);...
    cos(0) cos(2*a) cos(4*a) cos(6*a) cos(8*a) cos(10*a) cos(12*a);...
    sin(0) sin(2*a) sin(4*a) sin(6*a) sin(8*a) sin(10*a) sin(12*a);...
    cos(0) cos(3*a) cos(6*a) cos(9*a) cos(12*a) cos(15*a) cos(18*a);...
    sin(0) sin(3*a) sin(6*a) sin(9*a) sin(12*a) sin(15*a) sin(18*a);...
    1 1 1 1 1 1 1];
sindqxy=mapping*[s1;s2;s3;s4;s5;s6;s7];
sindq1x=sindqxy(1,:);
sindq1y=sindqxy(2,:);
sindq2x=sindqxy(3,:);
sindq2y=sindqxy(4,:);
sindq3x=sindqxy(5,:);
sindq3y=sindqxy(6,:);
angle=atan2(sindq1y,sindq1x);
for i=1:1:fsw*time
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
% figure;%陪ボsin计
% plot(sindq1y,sindq1x);
% hold on;
A = [6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7];
DQ=zeros(7,126);
M1=zeros(2,126);
M2=zeros(2,126);
M3=zeros(2,126);
for i=1:126
    binary=dec2bin(i,7);
    bin=binary-48;
    DQ(:,i)=mapping*A*[bin(1);bin(2);bin(3);bin(4);bin(5);bin(6);bin(7)];
    M1(:,i)=[DQ(1,i);DQ(2,i)];
    M2(:,i)=[DQ(3,i);DQ(4,i)];
    M3(:,i)=[DQ(5,i);DQ(6,i)];
end
% quiver(zeros(1,126),zeros(1,126),DQ(1,:),DQ(2,:))
% title('dq-plane 1'); 
% for i=1:126
%     k=num2str(i);
%     text(DQ(1,i),DQ(2,i),k);
% end
% figure;%陪ボ计
% quiver(zeros(1,126),zeros(1,126),DQ(3,:),DQ(4,:))
% title('dq-plane 2'); 
% for i=1:126
%     k=num2str(i);
%     text(DQ(3,i),DQ(4,i),k);
% end
% figure;%陪ボ计
% quiver(zeros(1,126),zeros(1,126),DQ(5,:),DQ(6,:))
% title('dq-plane 3'); 
% for i=1:126
%     k=num2str(i);
%     text(DQ(5,i),DQ(6,i),k);
% end

duty=zeros(6,fsw*time);
sector=zeros(1,fsw*time);
v1=zeros(7,fsw*time);
v2=zeros(7,fsw*time);
v3=zeros(7,fsw*time);
v4=zeros(7,fsw*time);
v5=zeros(7,fsw*time);
v6=zeros(7,fsw*time);
for i=1:fsw*time%耞sector
if angle(i)>0 && angle(i)<=pi/7;%sector1
    sector(i)=1;
    duty(:,i)=[M1(:,64) M1(:,96) M1(:,97) M1(:,113) M1(:,115) M1(:,123);...
        M2(:,64) M2(:,96) M2(:,97) M2(:,113) M2(:,115) M2(:,123);...
        M3(:,64) M3(:,96) M3(:,97) M3(:,113) M3(:,115) M3(:,123)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[1;0;0;0;0;0;0];
    v2(:,i)=[1;1;0;0;0;0;0];
    v3(:,i)=[1;1;0;0;0;0;1];
    v4(:,i)=[1;1;1;0;0;0;1];
    v5(:,i)=[1;1;1;0;0;1;1];
    v6(:,i)=[1;1;1;1;0;1;1];
end
if angle(i)>pi/7 && angle(i)<=2*pi/7;%sector2
    sector(i)=2;
    duty(:,i)=[M1(:,32) M1(:,96) M1(:,112) M1(:,113) M1(:,121) M1(:,123);...
        M2(:,32) M2(:,96) M2(:,112) M2(:,113) M2(:,121) M2(:,123);...
        M3(:,32) M3(:,96) M3(:,112) M3(:,113) M3(:,121) M3(:,123)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;1;0;0;0;0;0];
    v2(:,i)=[1;1;0;0;0;0;0];
    v3(:,i)=[1;1;1;0;0;0;0];
    v4(:,i)=[1;1;1;0;0;0;1];
    v5(:,i)=[1;1;1;1;0;0;1];
    v6(:,i)=[1;1;1;1;0;1;1];
end   
if angle(i)>2*pi/7 && angle(i)<=3*pi/7;%sector3
    sector(i)=3;
    duty(:,i)=[M1(:,32) M1(:,48) M1(:,112) M1(:,120) M1(:,121) M1(:,125);...
        M2(:,32) M2(:,48) M2(:,112) M2(:,120) M2(:,121) M2(:,125);...
        M3(:,32) M3(:,48) M3(:,112) M3(:,120) M3(:,121) M3(:,125)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;1;0;0;0;0;0];
    v2(:,i)=[0;1;1;0;0;0;0];
    v3(:,i)=[1;1;1;0;0;0;0];
    v4(:,i)=[1;1;1;1;0;0;0];
    v5(:,i)=[1;1;1;1;0;0;1];
    v6(:,i)=[1;1;1;1;1;0;1];
end
if angle(i)>3*pi/7 && angle(i)<=4*pi/7;%sector4
    sector(i)=4;
    duty(:,i)=[M1(:,16) M1(:,48) M1(:,56) M1(:,120) M1(:,124) M1(:,125);...
        M2(:,16) M2(:,48) M2(:,56) M2(:,120) M2(:,124) M2(:,125);...
        M3(:,16) M3(:,48) M3(:,56) M3(:,120) M3(:,124) M3(:,125)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;1;0;0;0;0];
    v2(:,i)=[0;1;1;0;0;0;0];
    v3(:,i)=[0;1;1;1;0;0;0];
    v4(:,i)=[1;1;1;1;0;0;0];
    v5(:,i)=[1;1;1;1;1;0;0];
    v6(:,i)=[1;1;1;1;1;0;1];
end
if angle(i)>4*pi/7 && angle(i)<=5*pi/7;%sector5
    sector(i)=5;
    duty(:,i)=[M1(:,16) M1(:,24) M1(:,56) M1(:,60) M1(:,124) M1(:,126);...
        M2(:,16) M2(:,24) M2(:,56) M2(:,60) M2(:,124) M2(:,126);...
        M3(:,16) M3(:,24) M3(:,56) M3(:,60) M3(:,124) M3(:,126)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;1;0;0;0;0];
    v2(:,i)=[0;0;1;1;0;0;0];
    v3(:,i)=[0;1;1;1;0;0;0];
    v4(:,i)=[0;1;1;1;1;0;0];
    v5(:,i)=[1;1;1;1;1;0;0];
    v6(:,i)=[1;1;1;1;1;1;0];
end 
if angle(i)>5*pi/7 && angle(i)<=6*pi/7;%sector6
    sector(i)=6;
    duty(:,i)=[M1(:,8) M1(:,24) M1(:,28) M1(:,60) M1(:,62) M1(:,126);...
        M2(:,8) M2(:,24) M2(:,28) M2(:,60) M2(:,62) M2(:,126);...
        M3(:,8) M3(:,24) M3(:,28) M3(:,60) M3(:,62) M3(:,126)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;1;0;0;0];
    v2(:,i)=[0;0;1;1;0;0;0];
    v3(:,i)=[0;0;1;1;1;0;0];
    v4(:,i)=[0;1;1;1;1;0;0];
    v5(:,i)=[0;1;1;1;1;1;0];
    v6(:,i)=[1;1;1;1;1;1;0];
end
if angle(i)>6*pi/7 && angle(i)<=pi;%sector7
    sector(i)=7;
    duty(:,i)=[M1(:,8) M1(:,12) M1(:,28) M1(:,30) M1(:,62) M1(:,63);...
        M2(:,8) M2(:,12) M2(:,28) M2(:,30) M2(:,62) M2(:,63);...
        M3(:,8) M3(:,12) M3(:,28) M3(:,30) M3(:,62) M3(:,63)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;1;0;0;0];
    v2(:,i)=[0;0;0;1;1;0;0];
    v3(:,i)=[0;0;1;1;1;0;0];
    v4(:,i)=[0;0;1;1;1;1;0];
    v5(:,i)=[0;1;1;1;1;1;0];
    v6(:,i)=[0;1;1;1;1;1;1];
end
if angle(i)>pi && angle(i)<=8*pi/7;%sector8
    sector(i)=8;
    duty(:,i)=[M1(:,4) M1(:,12) M1(:,14) M1(:,30) M1(:,31) M1(:,63);...
        M2(:,4) M2(:,12) M2(:,14) M2(:,30) M2(:,31) M2(:,63);...
        M3(:,4) M3(:,12) M3(:,14) M3(:,30) M3(:,31) M3(:,63)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;0;1;0;0];
    v2(:,i)=[0;0;0;1;1;0;0];
    v3(:,i)=[0;0;0;1;1;1;0];
    v4(:,i)=[0;0;1;1;1;1;0];
    v5(:,i)=[0;0;1;1;1;1;1];
    v6(:,i)=[0;1;1;1;1;1;1];
end
if angle(i)>8*pi/7 && angle(i)<=9*pi/7;%sector9
    sector(i)=9;
    duty(:,i)=[M1(:,4) M1(:,6) M1(:,14) M1(:,15) M1(:,31) M1(:,95);...
        M2(:,4) M2(:,6) M2(:,14) M2(:,15) M2(:,31) M2(:,95);...
        M3(:,4) M3(:,6) M3(:,14) M3(:,15) M3(:,31) M3(:,95)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;0;1;0;0];
    v2(:,i)=[0;0;0;0;1;1;0];
    v3(:,i)=[0;0;0;1;1;1;0];
    v4(:,i)=[0;0;0;1;1;1;1];
    v5(:,i)=[0;0;1;1;1;1;1];
    v6(:,i)=[1;0;1;1;1;1;1];
end
if angle(i)>9*pi/7 && angle(i)<=10*pi/7;%sector10
    sector(i)=10;
    duty(:,i)=[M1(:,2) M1(:,6) M1(:,7) M1(:,15) M1(:,79) M1(:,95);...
        M2(:,2) M2(:,6) M2(:,7) M2(:,15) M2(:,79) M2(:,95);...
        M3(:,2) M3(:,6) M3(:,7) M3(:,15) M3(:,79) M3(:,95)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;0;0;1;0];
    v2(:,i)=[0;0;0;0;1;1;0];
    v3(:,i)=[0;0;0;0;1;1;1];
    v4(:,i)=[0;0;0;1;1;1;1];
    v5(:,i)=[1;0;0;1;1;1;1];
    v6(:,i)=[1;0;1;1;1;1;1];
end
if angle(i)>10*pi/7 && angle(i)<=11*pi/7;%sector11
    sector(i)=11;
    duty(:,i)=[M1(:,2) M1(:,3) M1(:,7) M1(:,71) M1(:,79) M1(:,111);...
        M2(:,2) M2(:,3) M2(:,7) M2(:,71) M2(:,79) M2(:,111);...
        M3(:,2) M3(:,3) M3(:,7) M3(:,71) M3(:,79) M3(:,111)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;0;0;1;0];
    v2(:,i)=[0;0;0;0;0;1;1];
    v3(:,i)=[0;0;0;0;1;1;1];
    v4(:,i)=[1;0;0;0;1;1;1];
    v5(:,i)=[1;0;0;1;1;1;1];
    v6(:,i)=[1;1;0;1;1;1;1];
end
if angle(i)>11*pi/7 && angle(i)<=12*pi/7;%sector12
    sector(i)=12;
    duty(:,i)=[M1(:,1) M1(:,3) M1(:,67) M1(:,71) M1(:,103) M1(:,111);...
        M2(:,1) M2(:,3) M2(:,67) M2(:,71) M2(:,103) M2(:,111);...
        M3(:,1) M3(:,3) M3(:,67) M3(:,71) M3(:,103) M3(:,111)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;0;0;0;1];
    v2(:,i)=[0;0;0;0;0;1;1];
    v3(:,i)=[1;0;0;0;0;1;1];
    v4(:,i)=[1;0;0;0;1;1;1];
    v5(:,i)=[1;1;0;0;1;1;1];
    v6(:,i)=[1;1;0;1;1;1;1];
end
if angle(i)>12*pi/7 && angle(i)<=13*pi/7;%sector13
    sector(i)=13;
    duty(:,i)=[M1(:,1) M1(:,65) M1(:,67) M1(:,99) M1(:,103) M1(:,119);...
        M2(:,1) M2(:,65) M2(:,67) M2(:,99) M2(:,103) M2(:,119);...
        M3(:,1) M3(:,65) M3(:,67) M3(:,99) M3(:,103) M3(:,119)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[0;0;0;0;0;0;1];
    v2(:,i)=[1;0;0;0;0;0;1];
    v3(:,i)=[1;0;0;0;0;1;1];
    v4(:,i)=[1;1;0;0;0;1;1];
    v5(:,i)=[1;1;0;0;1;1;1];
    v6(:,i)=[1;1;1;0;1;1;1];
end
if angle(i)>13*pi/7 && angle(i)<=2*pi;%sector14
    sector(i)=14;
    duty(:,i)=[M1(:,64) M1(:,65) M1(:,97) M1(:,99) M1(:,115) M1(:,119);...
        M2(:,64) M2(:,65) M2(:,97) M2(:,99) M2(:,115) M2(:,119);...
        M3(:,64) M3(:,65) M3(:,97) M3(:,99) M3(:,115) M3(:,119)]\...
        [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    v1(:,i)=[1;0;0;0;0;0;0];
    v2(:,i)=[1;0;0;0;0;0;1];
    v3(:,i)=[1;1;0;0;0;0;1];
    v4(:,i)=[1;1;0;0;0;1;1];
    v5(:,i)=[1;1;1;0;0;1;1];
    v6(:,i)=[1;1;1;0;1;1;1];
end
end
a1=duty(1,:);
a2=duty(2,:);
a3=duty(3,:);
a4=duty(4,:);
a5=duty(5,:);
a6=duty(6,:);
az=1-a1-a2-a3-a4-a5-a6;
a0=az/2;
a7=a0;
% k=[sector;a0;a1;a2;a3;a4;a5;a6;a7];

T1=a1;%礚&Τ弘非
T2=a2;
T3=a3;
T4=a4;
T5=a5;
T6=a6;
Tz=az;
Q1=zeros(1,fsw*time);
Q2=zeros(1,fsw*time);
Q3=zeros(1,fsw*time);
Q4=zeros(1,fsw*time);
Q5=zeros(1,fsw*time);
Q6=zeros(1,fsw*time);
Qz=zeros(1,fsw*time);
D1=zeros(1,fsw*time);
D2=zeros(1,fsw*time);
D3=zeros(1,fsw*time);
D4=zeros(1,fsw*time);
D5=zeros(1,fsw*time);
D6=zeros(1,fsw*time);
Dz=zeros(1,fsw*time);
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(i),2)==1%1,sector碞琌计0,sector碞琌案计
        Q1(i)=((M1(1,64))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T1(i);
        Q2(i)=((M1(1,115))*cos((pi/7)*sector(i)-angle(i))-amp)*T2(i);
        Q3(i)=((M1(1,97))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T3(i);
        Q4(i)=((M1(1,97))*cos((pi/7)*sector(i)-angle(i))-amp)*T4(i);
        Q5(i)=((M1(1,115))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T5(i);
        Q6(i)=((M1(1,64))*cos((pi/7)*sector(i)-angle(i))-amp)*T6(i);
        Qz(i)=-amp*Tz(i);
        D1(i)=((M1(1,64))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T1(i);
        D2(i)=((-M1(1,115))*sin((pi/7)*sector(i)-angle(i)))*T2(i);
        D3(i)=((M1(1,97))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T3(i);
        D4(i)=((-M1(1,97))*sin((pi/7)*sector(i)-angle(i)))*T4(i);
        D5(i)=((M1(1,115))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T5(i);
        D6(i)=((-M1(1,64))*sin((pi/7)*sector(i)-angle(i)))*T6(i);
        Dz(i)=0;
    else
        Q1(i)=((M1(1,64))*cos((pi/7)*sector(i)-angle(i))-amp)*T1(i);
        Q2(i)=((M1(1,115))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T2(i);
        Q3(i)=((M1(1,97))*cos((pi/7)*sector(i)-angle(i))-amp)*T3(i);
        Q4(i)=((M1(1,97))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T4(i);
        Q5(i)=((M1(1,115))*cos((pi/7)*sector(i)-angle(i))-amp)*T5(i);
        Q6(i)=((M1(1,64))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T6(i);
        Qz(i)=-amp*Tz(i); 
        D1(i)=((M1(1,64))*sin((pi/7)*sector(i)-angle(i)))*T1(i);
        D2(i)=((-M1(1,115))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T2(i);
        D3(i)=((M1(1,97))*sin((pi/7)*sector(i)-angle(i)))*T3(i);
        D4(i)=((-M1(1,97))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T4(i);
        D5(i)=((M1(1,115))*sin((pi/7)*sector(i)-angle(i)))*T5(i);
        D6(i)=((-M1(1,64))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T6(i);
        Dz(i)=0;
    end
end
QQQQ=Qz+Q1+Q2+Q3+Q4+Q5+Q6;
DDDD=Dz+D1+D2+D3+D4+D5+D6;
% figure;
% plot(angle*180/pi,Q1,'r');
% hold on;
% plot(angle*180/pi,Q2,'b');
% hold on;
% plot(angle*180/pi,Q3,'g');
% hold on;
% plot(angle*180/pi,Q4,'y');
% hold on;
% plot(angle*180/pi,Q5,'k');
% hold on;
% plot(angle*180/pi,Q6,'c');
%{
HDFQ1=zeros(1,fsw*time);
HDFQ2=zeros(1,fsw*time);
HDFQ3=zeros(1,fsw*time);
HDFQ4=zeros(1,fsw*time);
HDFQ5=zeros(1,fsw*time);
HDFQ6=zeros(1,fsw*time);
HDFQ7=zeros(1,fsw*time);
HDFQ8=zeros(1,fsw*time);
HDFD1=zeros(1,fsw*time);
HDFD2=zeros(1,fsw*time);
HDFD3=zeros(1,fsw*time);
HDFD4=zeros(1,fsw*time);
HDFD5=zeros(1,fsw*time);
HDFD6=zeros(1,fsw*time);
HDFD7=zeros(1,fsw*time);
HDFD8=zeros(1,fsw*time);
for i=1:1:fsw*time%HDF01234567
    P1=0;
    P2=Qz(i)/2;
    P3=Qz(i)/2+Q1(i);
    P4=Qz(i)/2+Q1(i)+Q2(i);
    P5=Qz(i)/2+Q1(i)+Q2(i)+Q3(i);
    P6=Qz(i)/2+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    P7=Qz(i)/2+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i);
    P8=-Qz(i)/2;
    HDFQ1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2;
    HDFQ2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDFQ3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDFQ4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T3(i);
    HDFQ5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T4(i);
    HDFQ6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T5(i);
    HDFQ7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T6(i);
    HDFQ8(i)=((P8)^2+(P8)*(P1)+(P1)^2)*Tz(i)/2;
% HDFD1(i)=((0)^2+(0)*(Dz(i)/2)+(Dz(i)/2)^2)*Tz(i)/2;
% HDFD2(i)=((Dz(i)/2)^2+(Dz(i)/2)*(Dz(i)/2+D1(i))+(Dz(i)/2+D1(i))^2)*T1(i);
% HDFD3(i)=((Dz(i)/2+D1(i))^2+(Dz(i)/2+D1(i))*(Dz(i)/2+D1(i)+D2(i))+(Dz(i)/2+D1(i)+D2(i))^2)*T2(i);
% HDFD4(i)=((Dz(i)/2+D1(i)+D2(i))^2+(Dz(i)/2+D1(i)+D2(i))*(Dz(i)/2+D1(i)+D2(i)+D3(i))+(Dz(i)/2+D1(i)+D2(i)+D3(i))^2)*T3(i);
% HDFD5(i)=((Dz(i)/2+D1(i)+D2(i)+D3(i))^2+(Dz(i)/2+D1(i)+D2(i)+D3(i))*(Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i))+(Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i))^2)*T4(i);
% HDFD6(i)=((Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i))^2+(Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i))*(Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i)+D5(i))+(Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i)+D5(i))^2)*T5(i);
% HDFD7(i)=((Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i)+D5(i))^2+(Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i)+D5(i))*(-Dz(i)/2)+(-Dz(i)/2)^2)*T6(i);
% HDFD8(i)=((-Dz(i)/2)^2+(-Dz(i)/2)*(0)+(0)^2)*Tz(i)/2;
end
HDFQ=HDFQ1+HDFQ2+HDFQ3+HDFQ4+HDFQ5+HDFQ6+HDFQ7+HDFQ8;
HDFD=HDFD1+HDFD2+HDFD3+HDFD4+HDFD5+HDFD6+HDFD7+HDFD8;
HDF=HDFQ+HDFD;
HDF0121Q1=zeros(1,fsw*time);
HDF0121Q2=zeros(1,fsw*time);
HDF0121Q3=zeros(1,fsw*time);
HDF0121Q4=zeros(1,fsw*time);
HDF0121Q5=zeros(1,fsw*time);
HDF0121Q6=zeros(1,fsw*time);
HDF0121Q7=zeros(1,fsw*time);
HDF0121Q8=zeros(1,fsw*time);
HDF0121D1=zeros(1,fsw*time);
HDF0121D2=zeros(1,fsw*time);
HDF0121D3=zeros(1,fsw*time);
HDF0121D4=zeros(1,fsw*time);
HDF0121D5=zeros(1,fsw*time);
HDF0121D6=zeros(1,fsw*time);
HDF0121D7=zeros(1,fsw*time);
HDF0121D8=zeros(1,fsw*time);
for i=1:1:fsw*time%HDF01234565
    P1=0;
    P2=Qz(i);
    P3=Qz(i)+Q1(i);
    P4=Qz(i)+Q1(i)+Q2(i);
    P5=Qz(i)+Q1(i)+Q2(i)+Q3(i);
    P6=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    P7=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i)*fraction3;
    P8=-Q5(i)*(1-fraction3);
    HDF0121Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDF0121Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDF0121Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T3(i);
    HDF0121Q5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T4(i);
    HDF0121Q6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3;
    HDF0121Q7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T6(i);
    HDF0121Q8(i)=((P8)^2+(P8)*(P1)+(P1)^2)*T5(i)*(1-fraction3);
% HDF0121D1(i)=((0)^2+(0)*(Dz(i))+(Dz(i))^2)*Tz(i);
% HDF0121D2(i)=((Dz(i))^2+(Dz(i))*(Dz(i)+D1(i))+(Dz(i)+D1(i))^2)*T1(i);
% HDF0121D3(i)=((Dz(i)+D1(i))^2+(Dz(i)+D1(i))*(Dz(i)+D1(i)+D2(i))+(Dz(i)+D1(i)+D2(i))^2)*T2(i);
% HDF0121D4(i)=((Dz(i)+D1(i)+D2(i))^2+(Dz(i)+D1(i)+D2(i))*(Dz(i)+D1(i)+D2(i)+D3(i))+(Dz(i)+D1(i)+D2(i)+D3(i))^2)*T3(i);
% HDF0121D5(i)=((Dz(i)+D1(i)+D2(i)+D3(i))^2+(Dz(i)+D1(i)+D2(i)+D3(i))*(Dz(i)+D1(i)+D2(i)+D3(i)+D4(i))+(Dz(i)+D1(i)+D2(i)+D3(i)+D4(i))^2)*T4(i);
% HDF0121D6(i)=((Dz(i)+D1(i)+D2(i)+D3(i)+D4(i))^2+(Dz(i)+D1(i)+D2(i)+D3(i)+D4(i))*(Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)*fraction3)+(Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)*fraction3)^2)*T5(i)*fraction3;
% HDF0121D7(i)=((Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)*fraction3)^2+(Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)*fraction3)*(-D5(i)*(1-fraction3))+(-D5(i)*(1-fraction3))^2)*T6(i);
% HDF0121D8(i)=((-D5(i)*(1-fraction3))^2+(-D5(i)*(1-fraction3))*(0)+(0)^2)*T5(i)*(1-fraction3);
end
HDF0121Q=HDF0121Q1+HDF0121Q2+HDF0121Q3+HDF0121Q4+HDF0121Q5+HDF0121Q6+HDF0121Q7+HDF0121Q8;
HDF0121D=HDF0121D1+HDF0121D2+HDF0121D3+HDF0121D4+HDF0121D5+HDF0121D6+HDF0121D7+HDF0121D8;
HDF0121=HDF0121Q+HDF0121D;
HDF7212Q1=zeros(1,fsw*time);
HDF7212Q2=zeros(1,fsw*time);
HDF7212Q3=zeros(1,fsw*time);
HDF7212Q4=zeros(1,fsw*time);
HDF7212Q5=zeros(1,fsw*time);
HDF7212Q6=zeros(1,fsw*time);
HDF7212Q7=zeros(1,fsw*time);
HDF7212Q8=zeros(1,fsw*time);
HDF7212D1=zeros(1,fsw*time);
HDF7212D2=zeros(1,fsw*time);
HDF7212D3=zeros(1,fsw*time);
HDF7212D4=zeros(1,fsw*time);
HDF7212D5=zeros(1,fsw*time);
HDF7212D6=zeros(1,fsw*time);
HDF7212D7=zeros(1,fsw*time);
HDF7212D8=zeros(1,fsw*time);
for i=1:1:fsw*time%HDF76543212
    P1=0;
    P2=Qz(i);
    P3=Qz(i)+Q6(i);
    P4=Qz(i)+Q6(i)+Q5(i);
    P5=Qz(i)+Q6(i)+Q5(i)+Q4(i);
    P6=Qz(i)+Q6(i)+Q5(i)+Q4(i)+Q3(i);
    P7=Qz(i)+Q6(i)+Q5(i)+Q4(i)+Q3(i)+Q2(i)*fraction3;
    P8=-Q2(i)*(1-fraction3);
    HDF0121Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T6(i);
    HDF0121Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T5(i);
    HDF0121Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T4(i);
    HDF0121Q5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T3(i);
    HDF0121Q6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3;
    HDF0121Q7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T1(i);
    HDF0121Q8(i)=((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction3);
% HDF7212D1(i)=((0)^2+(0)*(Dz(i))+(Dz(i))^2)*Tz(i);
% HDF7212D2(i)=((Dz(i))^2+(Dz(i))*(Dz(i)+D6(i))+(Dz(i)+D6(i))^2)*T6(i);
% HDF7212D3(i)=((Dz(i)+D6(i))^2+(Dz(i)+D6(i))*(Dz(i)+D6(i)+D5(i))+(Dz(i)+D6(i)+D5(i))^2)*T5(i);
% HDF7212D4(i)=((Dz(i)+D6(i)+D5(i))^2+(Dz(i)+D6(i)+D5(i))*(Dz(i)+D6(i)+D5(i)+D4(i))+(Dz(i)+D6(i)+D5(i)+D4(i))^2)*T4(i);
% HDF7212D5(i)=((Dz(i)+D6(i)+D5(i)+D4(i))^2+(Dz(i)+D6(i)+D5(i)+D4(i))*(Dz(i)+D6(i)+D5(i)+D4(i)+D3(i))+(Dz(i)+D6(i)+D5(i)+D4(i)+D3(i))^2)*T3(i);
% HDF7212D6(i)=((Dz(i)+D6(i)+D5(i)+D4(i)+D3(i))^2+(Dz(i)+D6(i)+D5(i)+D4(i)+D3(i))*(Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)*fraction3)+(Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)*fraction3)^2)*T2(i)*fraction3;
% HDF7212D7(i)=((Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)*fraction3)^2+(Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)*fraction3)*(-D2(i)*(1-fraction3))+(-D2(i)*(1-fraction3))^2)*T1(i);
% HDF7212D8(i)=((-D2(i)*(1-fraction3))^2+(-D2(i)*(1-fraction3))*(0)+(0)^2)*T2(i)*(1-fraction3);
end
HDF7212Q=HDF7212Q1+HDF7212Q2+HDF7212Q3+HDF7212Q4+HDF7212Q5+HDF7212Q6+HDF7212Q7+HDF7212Q8;
HDF7212D=HDF7212D1+HDF7212D2+HDF7212D3+HDF7212D4+HDF7212D5+HDF7212D6+HDF7212D7+HDF7212D8;
HDF7212=HDF7212Q+HDF7212D;
%}
for i=1:1:fsw*time%01234567
    P1=0;
    P2=Qz(i)/2;
    P3=Qz(i)/2+Q1(i);
    P4=Qz(i)/2+Q1(i)+Q2(i);
    P5=Qz(i)/2+Q1(i)+Q2(i)+Q3(i);
    P6=Qz(i)/2+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    P7=Qz(i)/2+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i);
    P8=Qz(i)/2+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i)+Q6(i);
    P9=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i)+Q6(i);
    HDFQ1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2;
    HDFQ2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDFQ3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDFQ4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T3(i);
    HDFQ5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T4(i);
    HDFQ6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T5(i);
    HDFQ7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T6(i);
    HDFQ8(i)=((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Dz(i)/2;
    R3=Dz(i)/2+D1(i);
    R4=Dz(i)/2+D1(i)+D2(i);
    R5=Dz(i)/2+D1(i)+D2(i)+D3(i);
    R6=Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i);
    R7=Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i)+D5(i);
    R8=Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)+D6(i);
    R9=Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)+D6(i);
    HDFD1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2;
    HDFD2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i);
    HDFD3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDFD4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T3(i);
    HDFD5(i)=((R5)^2+(R5)*(R6)+(R6)^2)*T4(i);
    HDFD6(i)=((R6)^2+(R6)*(R7)+(R7)^2)*T5(i);
    HDFD7(i)=((R7)^2+(R7)*(R8)+(R8)^2)*T6(i);
    HDFD8(i)=((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Qz(i);
    P3=Qz(i)+Q1(i);
    P4=Qz(i)+Q1(i)+Q2(i);
    P5=Qz(i)+Q1(i)+Q2(i)+Q3(i);
    P6=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    P7=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i)*fraction3;
    P8=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i)*fraction3+Q6(i);
    P9=Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i)+Q5(i)+Q6(i);
    HDF0121Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDF0121Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDF0121Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T3(i);
    HDF0121Q5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T4(i);
    HDF0121Q6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3;
    HDF0121Q7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T6(i);
    HDF0121Q8(i)=((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Dz(i);
    R3=Dz(i)+D1(i);
    R4=Dz(i)+D1(i)+D2(i);
    R5=Dz(i)+D1(i)+D2(i)+D3(i);
    R6=Dz(i)+D1(i)+D2(i)+D3(i)+D4(i);
    R7=Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)*fraction3;
    R8=Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)*fraction3+D6(i);
    R9=Dz(i)+D1(i)+D2(i)+D3(i)+D4(i)+D5(i)+D6(i);
    HDF0121D1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF0121D2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i);
    HDF0121D3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDF0121D4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T3(i);
    HDF0121D5(i)=((R5)^2+(R5)*(R6)+(R6)^2)*T4(i);
    HDF0121D6(i)=((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3;
    HDF0121D7(i)=((R7)^2+(R7)*(R8)+(R8)^2)*T6(i);
    HDF0121D8(i)=((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz(i);
    P3=Qz(i)+Q6(i);
    P4=Qz(i)+Q6(i)+Q5(i);
    P5=Qz(i)+Q6(i)+Q5(i)+Q4(i);
    P6=Qz(i)+Q6(i)+Q5(i)+Q4(i)+Q3(i);
    P7=Qz(i)+Q6(i)+Q5(i)+Q4(i)+Q3(i)+Q2(i)*fraction3;
    P8=Qz(i)+Q6(i)+Q5(i)+Q4(i)+Q3(i)+Q2(i)*fraction3+Q1(i);
    P9=Qz(i)+Q6(i)+Q5(i)+Q4(i)+Q3(i)+Q2(i)+Q1(i);
    HDF7212Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF7212Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T6(i);
    HDF7212Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T5(i);
    HDF7212Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T4(i);
    HDF7212Q5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T3(i);
    HDF7212Q6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3;
    HDF7212Q7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T1(i);
    HDF7212Q8(i)=((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Dz(i);
    R3=Dz(i)+D6(i);
    R4=Dz(i)+D6(i)+D5(i);
    R5=Dz(i)+D6(i)+D5(i)+D4(i);
    R6=Dz(i)+D6(i)+D5(i)+D4(i)+D3(i);
    R7=Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)*fraction3;
    R8=Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)*fraction3+D1(i);
    R9=Dz(i)+D6(i)+D5(i)+D4(i)+D3(i)+D2(i)+D1(i);
    HDF7212D1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF7212D2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T6(i);
    HDF7212D3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T5(i);
    HDF7212D4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T4(i);
    HDF7212D5(i)=((R5)^2+(R5)*(R6)+(R6)^2)*T3(i);
    HDF7212D6(i)=((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3;
    HDF7212D7(i)=((R7)^2+(R7)*(R8)+(R8)^2)*T1(i);
    HDF7212D8(i)=((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction3);
end
HDFQ=HDFQ1+HDFQ2+HDFQ3+HDFQ4+HDFQ5+HDFQ6+HDFQ7+HDFQ8;
HDFD=HDFD1+HDFD2+HDFD3+HDFD4+HDFD5+HDFD6+HDFD7+HDFD8;
HDF=HDFQ+HDFD;
HDF0121Q=HDF0121Q1+HDF0121Q2+HDF0121Q3+HDF0121Q4+HDF0121Q5+HDF0121Q6+HDF0121Q7+HDF0121Q8;
HDF0121D=HDF0121D1+HDF0121D2+HDF0121D3+HDF0121D4+HDF0121D5+HDF0121D6+HDF0121D7+HDF0121D8;
HDF0121=HDF0121Q+HDF0121D;
HDF7212Q=HDF7212Q1+HDF7212Q2+HDF7212Q3+HDF7212Q4+HDF7212Q5+HDF7212Q6+HDF7212Q7+HDF7212Q8;
HDF7212D=HDF7212D1+HDF7212D2+HDF7212D3+HDF7212D4+HDF7212D5+HDF7212D6+HDF7212D7+HDF7212D8;
HDF7212=HDF7212Q+HDF7212D;

threeHDF=[HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
allarea3=sum(abs(HDF-minHDFvalue3));

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
T5bit=round(T5*bit/2);
T6bit=round(T6*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit-T5bit-T6bit;
TTT=T1bit+T2bit+T3bit+T4bit+T5bit+T6bit+Tzbit;
TTTTT=[T1bit;T2bit;T3bit;T4bit;Tzbit;TTT];
out=zeros(7,fsw*time);
for i=1:fsw*time%璶祘Α絏
   switch minHDFcase3(i)
       case 1%01234567
           out0=repmat([0;0;0;0;0;0;0],1,floor(Tzbit(i)/2));
           out1=repmat(v1(:,i),1,T1bit(i));
           out2=repmat(v2(:,i),1,T2bit(i));
           out3=repmat(v3(:,i),1,T3bit(i));
           out4=repmat(v4(:,i),1,T4bit(i));
           out5=repmat(v5(:,i),1,T5bit(i));
           out6=repmat(v6(:,i),1,T6bit(i));
           out7=repmat([1;1;1;1;1;1;1],1,Tzbit(i)-floor(Tzbit(i)/2));
           out01234567=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out01234567;
       case 2%01234565
           out0=repmat([0;0;0;0;0;0;0],1,Tzbit(i));
           out1=repmat(v1(:,i),1,T1bit(i));
           out2=repmat(v2(:,i),1,T2bit(i));
           out3=repmat(v3(:,i),1,T3bit(i));
           out4=repmat(v4(:,i),1,T4bit(i));
           out5=repmat(v5(:,i),1,floor(T5bit(i)*fraction3));
           out6=repmat(v6(:,i),1,T6bit(i));
           out7=repmat(v5(:,i),1,T5bit(i)-floor(T5bit(i)*fraction3));
           out01234565=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out01234565;
       case 3%76543212
           out0=repmat([1;1;1;1;1;1;1],1,Tzbit(i));
           out1=repmat(v6(:,i),1,T6bit(i));
           out2=repmat(v5(:,i),1,T5bit(i));
           out3=repmat(v4(:,i),1,T4bit(i));
           out4=repmat(v3(:,i),1,T3bit(i));
           out5=repmat(v2(:,i),1,floor(T2bit(i)*fraction3));
           out6=repmat(v1(:,i),1,T1bit(i));
           out7=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)*fraction3));
           out76543212=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out76543212;
   end
end
out_2(1,:)=out(1,2:end);
out_2(2,:)=out(2,2:end);
out_2(3,:)=out(3,2:end);
out_2(4,:)=out(4,2:end);
out_2(5,:)=out(5,2:end);
out_2(6,:)=out(6,2:end);
out_2(7,:)=out(7,2:end);
out_2(1,fsw*time*bit)=out(1,1);
out_2(2,fsw*time*bit)=out(2,1);
out_2(3,fsw*time*bit)=out(3,1);
out_2(4,fsw*time*bit)=out(4,1);
out_2(5,fsw*time*bit)=out(5,1);
out_2(6,fsw*time*bit)=out(6,1);
out_2(7,fsw*time*bit)=out(7,1);
out_3=[out;out_2];
diff1=sum(abs(out(1,:)-out_2(1,:)));
diff2=sum(abs(out(2,:)-out_2(2,:)));
diff3=sum(abs(out(3,:)-out_2(3,:)));
diff4=sum(abs(out(4,:)-out_2(4,:)));
diff5=sum(abs(out(5,:)-out_2(5,:)));
diff6=sum(abs(out(6,:)-out_2(6,:)));
diff7=sum(abs(out(7,:)-out_2(7,:)));
fprintf('out1ち传Ω计=%f\n',diff1);
% fprintf('out2ち传Ω计=%f\n',diff2);
% fprintf('out3ち传Ω计=%f\n',diff3);
% fprintf('out4ち传Ω计=%f\n',diff4);
% fprintf('out5ち传Ω计=%f\n',diff5);
% fprintf('out6ち传Ω计=%f\n',diff6);
% fprintf('out7ち传Ω计=%f\n',diff7);

%{
N = time*fsw*bit;
Van=[6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vbn=[-1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vcn=[-1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7]*out;
Vdn=[-1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7]*out;
Ven=[-1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7]*out;
Vfn=[-1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7]*out;
Vgn=[-1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7]*out;
fft_out1 = fft(Van,N);%6:0.499ㄢ狠常Τ
mag_out1 = abs(fft_out1)/(length(fft_out1)/2); 
fft_out2 = fft(Vbn,N); 
mag_out2 = abs(fft_out2)/(length(fft_out2)/2); 
fft_out3 = fft(Vcn,N); 
mag_out3 = abs(fft_out3)/(length(fft_out3)/2); 
fft_out4 = fft(Vdn,N);
mag_out4 = abs(fft_out4)/(length(fft_out4)/2); 
fft_out5 = fft(Ven,N);
mag_out5 = abs(fft_out5)/(length(fft_out5)/2); 
fft_out6 = fft(Vfn,N);
mag_out6 = abs(fft_out6)/(length(fft_out6)/2); 
fft_out7 = fft(Vgn,N);
mag_out7 = abs(fft_out7)/(length(fft_out7)/2); 
freq = (1:length(fft_out1)/2-1)*fsw*bit/length(fft_out1); 
%freq = (0:150000/2-1)*10;
draw1=mag_out1(2:1:N/2);%礶fft暗硂妓锣传
draw2=mag_out2(2:1:N/2);
draw3=mag_out3(2:1:N/2);
draw4=mag_out4(2:1:N/2);
draw5=mag_out5(2:1:N/2);
draw6=mag_out6(2:1:N/2);
draw7=mag_out7(2:1:N/2);
figure;
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
harm4 = mag_out4( 2*k+1: k: end/2 );
THD4 = 100*sqrt(  sum( harm4.^2)/mag_out4(k+1)^2 );
harm5 = mag_out5( 2*k+1: k: end/2 );
THD5 = 100*sqrt(  sum( harm5.^2)/mag_out5(k+1)^2 );
harm6 = mag_out6( 2*k+1: k: end/2 );
THD6 = 100*sqrt(  sum( harm6.^2)/mag_out6(k+1)^2 );
harm7 = mag_out7( 2*k+1: k: end/2 );
THD7 = 100*sqrt(  sum( harm7.^2)/mag_out7(k+1)^2 );
fprintf('VTHD1=%f\n',THD1);
fprintf('VTHD2=%f\n',THD2);
fprintf('VTHD3=%f\n',THD3);
fprintf('VTHD4=%f\n',THD4);
fprintf('VTHD5=%f\n',THD5);
fprintf('VTHD6=%f\n',THD6);
fprintf('VTHD7=%f\n',THD7);

i=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
mag_Van_1=mag_out1(1,i);
mag_Vbn_1=mag_out2(1,i);
mag_Vcn_1=mag_out3(1,i);
mag_Vdn_1=mag_out4(1,i);
mag_Ven_1=mag_out5(1,i);
mag_Vfn_1=mag_out6(1,i);
mag_Vgn_1=mag_out7(1,i);
figure;
semilogx(mag_Van_1,'r');
title('semilogx mag Van 1');
for i=1:fsw*bit*0.02/2-1
    mag_ithd_out1(i)=mag_Van_1(i)/(10*i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(10*i); 
    mag_ithd_out3(i)=mag_Vcn_1(i)/(10*i);  
    mag_ithd_out4(i)=mag_Vdn_1(i)/(10*i);  
    mag_ithd_out5(i)=mag_Ven_1(i)/(10*i);  
    mag_ithd_out6(i)=mag_Vfn_1(i)/(10*i);  
    mag_ithd_out7(i)=mag_Vgn_1(i)/(10*i);  
end
THD_Van=100*((sum(mag_ithd_out1.^2)-(mag_ithd_out1(1,1).^2))/(mag_ithd_out1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_ithd_out2.^2)-(mag_ithd_out2(1,1).^2))/(mag_ithd_out2(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_ithd_out3.^2)-(mag_ithd_out3(1,1).^2))/(mag_ithd_out3(1,1).^2)).^(1/2);
THD_Vdn=100*((sum(mag_ithd_out4.^2)-(mag_ithd_out4(1,1).^2))/(mag_ithd_out4(1,1).^2)).^(1/2);
THD_Ven=100*((sum(mag_ithd_out5.^2)-(mag_ithd_out5(1,1).^2))/(mag_ithd_out5(1,1).^2)).^(1/2);
THD_Vfn=100*((sum(mag_ithd_out6.^2)-(mag_ithd_out6(1,1).^2))/(mag_ithd_out6(1,1).^2)).^(1/2);
THD_Vgn=100*((sum(mag_ithd_out7.^2)-(mag_ithd_out7(1,1).^2))/(mag_ithd_out7(1,1).^2)).^(1/2);
fprintf('ITHD1=%f\n',THD_Van);
fprintf('ITHD2=%f\n',THD_Vbn);
fprintf('ITHD3=%f\n',THD_Vcn);
fprintf('ITHD4=%f\n',THD_Vdn);
fprintf('ITHD5=%f\n',THD_Ven);
fprintf('ITHD6=%f\n',THD_Vfn);
fprintf('ITHD7=%f\n',THD_Vgn);

A = [6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7];
B = [out(1,:);out(2,:);out(3,:);out(4,:);out(5,:);out(6,:);out(7,:)];  
C = 48*A*B; 
figure;
plot(C(1,1:1:fsw*bit*time));
hold on;
plot(C(2,1:1:fsw*bit*time));
hold on;
plot(C(3,1:1:fsw*bit*time));
hold on;
plot(C(4,1:1:fsw*bit*time));
hold on;
plot(C(5,1:1:fsw*bit*time));
hold on;
plot(C(6,1:1:fsw*bit*time));
hold on;
plot(C(7,1:1:fsw*bit*time));
title('筿溃');
% figure
% semilogx(C(1,1:1:fsw*bit*time))    
end_time=clock;
execution_time=end_time-start_time;
%}
N = time*fsw*bit;
Van=[6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vbn=[-1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vcn=[-1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7]*out;
Vdn=[-1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7]*out;
Ven=[-1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7]*out;
Vfn=[-1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7]*out;
Vgn=[-1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7]*out;
fft_Van=fft(Van,N)/N;%6:0.499ㄢ狠常Τ
mag_Van=abs(fft_Van)*2;
fft_Vbn=fft(Vbn,N)/N;
mag_Vbn=abs(fft_Vbn)*2;
fft_Vcn=fft(Vcn,N)/N;
mag_Vcn=abs(fft_Vcn)*2;
fft_Vdn=fft(Vdn,N)/N;
mag_Vdn=abs(fft_Vdn)*2;
fft_Ven=fft(Ven,N)/N;
mag_Ven=abs(fft_Ven)*2;
fft_Vfn=fft(Vfn,N)/N;
mag_Vfn=abs(fft_Vfn)*2;
fft_Vgn=fft(Vgn,N)/N;
mag_Vgn=abs(fft_Vgn)*2;
k=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:い丁
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
mag_Vfn_1=mag_Vfn(1,k);
mag_Vgn_1=mag_Vgn(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(繵^2/膀繵^2)
THD_Van=100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
THD_Vdn=100*((sum(mag_Vdn_1.^2)-(mag_Vdn_1(1,1).^2))/(mag_Vdn_1(1,1).^2)).^(1/2);
THD_Ven=100*((sum(mag_Ven_1.^2)-(mag_Ven_1(1,1).^2))/(mag_Ven_1(1,1).^2)).^(1/2);
THD_Vfn=100*((sum(mag_Vfn_1.^2)-(mag_Vfn_1(1,1).^2))/(mag_Vfn_1(1,1).^2)).^(1/2);
THD_Vgn=100*((sum(mag_Vgn_1.^2)-(mag_Vgn_1(1,1).^2))/(mag_Vgn_1(1,1).^2)).^(1/2);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n',THD_Vcn);
fprintf('VTHD4=%f\n',THD_Vdn);
fprintf('VTHD5=%f\n',THD_Ven);
fprintf('VTHD6=%f\n',THD_Vfn);
fprintf('VTHD7=%f\n',THD_Vgn);
%ITHD=繵瞯埃セō
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
mag_ithd_out6=zeros(1,N/5/2-1);
mag_ithd_out7=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%ソ搭埃そ畉+1=[(N/2-4)-6]/5+1
for i=1:length(mag_Van_1)
    mag_ithd_out1(i)=mag_Van_1(i)/(10*i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(10*i); 
    mag_ithd_out3(i)=mag_Vcn_1(i)/(10*i);  
    mag_ithd_out4(i)=mag_Vdn_1(i)/(10*i);  
    mag_ithd_out5(i)=mag_Ven_1(i)/(10*i);  
    mag_ithd_out6(i)=mag_Vfn_1(i)/(10*i);  
    mag_ithd_out7(i)=mag_Vgn_1(i)/(10*i);  
end
ITHD_Van=100*((sum(mag_ithd_out1.^2)-(mag_ithd_out1(1,1).^2))/(mag_ithd_out1(1,1).^2)).^(1/2);
ITHD_Vbn=100*((sum(mag_ithd_out2.^2)-(mag_ithd_out2(1,1).^2))/(mag_ithd_out2(1,1).^2)).^(1/2);
ITHD_Vcn=100*((sum(mag_ithd_out3.^2)-(mag_ithd_out3(1,1).^2))/(mag_ithd_out3(1,1).^2)).^(1/2);
ITHD_Vdn=100*((sum(mag_ithd_out4.^2)-(mag_ithd_out4(1,1).^2))/(mag_ithd_out4(1,1).^2)).^(1/2);
ITHD_Ven=100*((sum(mag_ithd_out5.^2)-(mag_ithd_out5(1,1).^2))/(mag_ithd_out5(1,1).^2)).^(1/2);
ITHD_Vfn=100*((sum(mag_ithd_out6.^2)-(mag_ithd_out6(1,1).^2))/(mag_ithd_out6(1,1).^2)).^(1/2);
ITHD_Vgn=100*((sum(mag_ithd_out7.^2)-(mag_ithd_out7(1,1).^2))/(mag_ithd_out7(1,1).^2)).^(1/2);
fprintf('ITHD1=%f\n',ITHD_Van);
fprintf('ITHD2=%f\n',ITHD_Vbn);
fprintf('ITHD3=%f\n',ITHD_Vcn);
fprintf('ITHD4=%f\n',ITHD_Vdn);
fprintf('ITHD5=%f\n',ITHD_Ven);
fprintf('ITHD6=%f\n',ITHD_Vfn);
fprintf('ITHD7=%f\n',ITHD_Vgn);
%赣衡,︽常衡ㄓ,钡ㄓ琌礶瓜
n2=2:1:N/2;%繷(繷材Τ粇箂)
mag_Van_2(1,n2-1)=mag_Van(1,n2);%5:0.499,10:0.00007,15...
mag_Vbn_2(1,n2-1)=mag_Vbn(1,n2); 
mag_Vcn_2(1,n2-1)=mag_Vcn(1,n2); 
mag_Vdn_2(1,n2-1)=mag_Vdn(1,n2); 
mag_Ven_2(1,n2-1)=mag_Ven(1,n2); 
mag_Vfn_2(1,n2-1)=mag_Vfn(1,n2); 
mag_Vgn_2(1,n2-1)=mag_Vgn(1,n2); 
mag_Van_3=zeros(1,(N/2-1)*10);
mag_Vbn_3=zeros(1,(N/2-1)*10);
mag_Vcn_3=zeros(1,(N/2-1)*10);
mag_Vdn_3=zeros(1,(N/2-1)*10);
mag_Ven_3=zeros(1,(N/2-1)*10);
mag_Vfn_3=zeros(1,(N/2-1)*10);
mag_Vgn_3=zeros(1,(N/2-1)*10);
for k=1:1:N/2-1;%5->50
mag_Van_3(1,(bit*fsw/N)*k)=mag_Van_2(1,k);
mag_Vbn_3(1,(bit*fsw/N)*k)=mag_Vbn_2(1,k);
mag_Vcn_3(1,(bit*fsw/N)*k)=mag_Vcn_2(1,k);
mag_Vdn_3(1,(bit*fsw/N)*k)=mag_Vdn_2(1,k);
mag_Ven_3(1,(bit*fsw/N)*k)=mag_Ven_2(1,k);
mag_Vfn_3(1,(bit*fsw/N)*k)=mag_Vfn_2(1,k);
mag_Vgn_3(1,(bit*fsw/N)*k)=mag_Vgn_2(1,k);
end
% figure===================================================================
% semilogx(mag_Van_3,'r'),title('mag Van 3');
end_time=clock;
execution_time=end_time-start_time;
ITHDavg=(ITHD_Van+ITHD_Vbn+ITHD_Vcn+ITHD_Vdn+ITHD_Ven+ITHD_Vfn+ITHD_Vgn)/7;
disperence=[point3;allarea3;THD_Van;ITHD_Van;ITHD_Vbn;ITHDavg];