clear%5/16ノ程HDF糶,蛾だガ
close all
clc
time = 0.02;
fsw = 15000;
bit = 400;
amp = 0.57;%0.57
fraction3 = 0.5;
f = 50;
t = 0:1/fsw:time-1/fsw;
interval = 0.01;
r = fix(amp/interval);
for j=1:1:r
   amp = interval*(j-1);
   s1(j,:) = amp*cos(2*pi*f*t);
   s2(j,:) = amp*cos(2*pi*f*t-2*pi/5);
   s3(j,:) = amp*cos(2*pi*f*t-4*pi/5);
   s4(j,:) = amp*cos(2*pi*f*t-6*pi/5);
   s5(j,:) = amp*cos(2*pi*f*t-8*pi/5);
end
a = 2*pi/5;
mapping=2/5*[...
    cos(0) cos(a) cos(2*a) cos(3*a) cos(4*a);...
    sin(0) sin(a) sin(2*a) sin(3*a) sin(4*a);...
    cos(0) cos(2*a) cos(4*a) cos(6*a) cos(8*a);...
    sin(0) sin(2*a) sin(4*a) sin(6*a) sin(8*a);...
    1 1 1 1 1];
for j=1:1:r
sindq(5*j-4,:) = mapping(1,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j-3,:) = mapping(2,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j-2,:) = mapping(3,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j-1,:) = mapping(4,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j,:) = mapping(5,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq_p1x(j,:) = sindq(5*j-4,:);
sindq_p1y(j,:) = sindq(5*j-3,:);
sindq_p2x(j,:) = sindq(5*j-2,:);
sindq_p2y(j,:) = sindq(5*j-1,:);
end
angle=atan2(sindq_p1y,sindq_p1x);
for i=1:1:fsw*time%==================================================================
    if angle(2,i)<0
        angle(2,i)=angle(2,i)+2*pi;
    end
end

for j=3:1:r
angle(j,:)=angle(2,:);
end
angle=angle+1e-10;
A = [4/5 -1/5 -1/5 -1/5 -1/5;...
    -1/5 4/5 -1/5 -1/5 -1/5;...
    -1/5 -1/5 4/5 -1/5 -1/5;...
    -1/5 -1/5 -1/5 4/5 -1/5;...
    -1/5 -1/5 -1/5 -1/5 4/5];
DQ=zeros(5,30);
M1=zeros(2,30);
M2=zeros(2,30);
for i=1:30
    binary=dec2bin(i,5);
    bin=binary-48;
    DQ(:,i)=mapping*A*[bin(1);bin(2);bin(3);bin(4);bin(5)];
    M1(:,i)=[DQ(1,i);DQ(2,i)];
    M2(:,i)=[DQ(3,i);DQ(4,i)];
end
for j=1:1:r
for i=1:fsw*time%耞sector
if angle(j,i)>0 && angle(j,i)<=pi/5;%sector1
    sector(j,i)=1;
    sb1=16;
    sb2=24;
    sb3=25;
    sb4=29;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>pi/5 && angle(j,i)<=2*pi/5;%sector2
    sector(j,i)=2;
    sb1=8;
    sb2=24;
    sb3=28;
    sb4=29;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end   
if angle(j,i)>2*pi/5 && angle(j,i)<=3*pi/5;%sector3
    sector(j,i)=3;
    sb1=8;
    sb2=12;
    sb3=28;
    sb4=30;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>3*pi/5 && angle(j,i)<=4*pi/5;%sector4
    sector(j,i)=4;
    sb1=4;
    sb2=12;
    sb3=14;
    sb4=30;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>4*pi/5 && angle(j,i)<=pi;%sector5
    sector(j,i)=5;
    sb1=4;
    sb2=6;
    sb3=14;
    sb4=15;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end 
if angle(j,i)>pi && angle(j,i)<=6*pi/5;%sector6
    sector(j,i)=6;
    sb1=2;
    sb2=6;
    sb3=7;
    sb4=15;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>6*pi/5 && angle(j,i)<=7*pi/5;%sector7
    sector(j,i)=7;
    sb1=2;
    sb2=3;
    sb3=7;
    sb4=23;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>7*pi/5 && angle(j,i)<=8*pi/5;%sector8
    sector(j,i)=8;
    sb1=1;
    sb2=3;
    sb3=19;
    sb4=23;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>8*pi/5 && angle(j,i)<=9*pi/5;%sector9
    sector(j,i)=9;
    sb1=1;
    sb2=17;
    sb3=19;
    sb4=27;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
if angle(j,i)>9*pi/5 && angle(j,i)<=2*pi;%sector10
    sector(j,i)=10;
    sb1=16;
    sb2=17;
    sb3=25;
    sb4=27;
    duty(:,i)=[M1(:,sb1) M1(:,sb2) M1(:,sb3) M1(:,sb4);...
              M2(:,sb1) M2(:,sb2) M2(:,sb3) M2(:,sb4)]\...
             [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i)];
end
Q5(j,i)=M2(1,sb1)*duty(1,i);
Q6(j,i)=M2(1,sb2)*duty(2,i);
Q7(j,i)=M2(1,sb3)*duty(3,i);
Q8(j,i)=M2(1,sb4)*duty(4,i);
Q9(j,i)=0;
D5(j,i)=M2(2,sb1)*duty(1,i);
D6(j,i)=M2(2,sb2)*duty(2,i);
D7(j,i)=M2(2,sb3)*duty(3,i);
D8(j,i)=M2(2,sb4)*duty(4,i);
D9(j,i)=0;
a1(j,i)=duty(1,i);
a2(j,i)=duty(2,i);
a3(j,i)=duty(3,i);
a4(j,i)=duty(4,i);
end
end
az=1-a1-a2-a3-a4;
a0=az/2;
a7=a0;

T1=a1;%礚蛤Τ弘非
T2=a2;
T3=a3;
T4=a4;
Tz=az;
for j=1:1:r
    amp=interval*(j-1);
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(j,i),2)==1%1,sector碞琌计0,sector碞琌案计
        Q1(j,i)=((M1(1,16))*cos(angle(j,i)-(pi/5)*(sector(j,i)-1))-amp)*T1(j,i);
        Q2(j,i)=((M1(1,25))*cos((pi/5)*sector(j,i)-angle(j,i))-amp)*T2(j,i);
        Q3(j,i)=((M1(1,25))*cos(angle(j,i)-(pi/5)*(sector(j,i)-1))-amp)*T3(j,i);
        Q4(j,i)=((M1(1,16))*cos((pi/5)*sector(j,i)-angle(j,i))-amp)*T4(j,i);
        Qz(j,i)=-amp*Tz(j,i);
        D1(j,i)=((M1(1,16))*sin(angle(j,i)-(pi/5)*(sector(j,i)-1)))*T1(j,i);
        D2(j,i)=((-M1(1,25))*sin((pi/5)*sector(j,i)-angle(j,i)))*T2(j,i);
        D3(j,i)=((M1(1,25))*sin(angle(j,i)-(pi/5)*(sector(j,i)-1)))*T3(j,i);
        D4(j,i)=((-M1(1,16))*sin((pi/5)*sector(j,i)-angle(j,i)))*T4(j,i);
        Dz(j,i)=0;
    else
        Q1(j,i)=((M1(1,16))*cos((pi/5)*sector(j,i)-angle(j,i))-amp)*T1(j,i);
        Q2(j,i)=((M1(1,25))*cos(angle(j,i)-(pi/5)*(sector(j,i)-1))-amp)*T2(j,i);
        Q3(j,i)=((M1(1,25))*cos((pi/5)*sector(j,i)-angle(j,i))-amp)*T3(j,i);
        Q4(j,i)=((M1(1,16))*cos(angle(j,i)-(pi/5)*(sector(j,i)-1))-amp)*T4(j,i);
        Qz(j,i)=-amp*Tz(j,i); 
        D1(j,i)=((M1(1,16))*sin((pi/5)*sector(j,i)-angle(j,i)))*T1(j,i);
        D2(j,i)=((-M1(1,25))*sin(angle(j,i)-(pi/5)*(sector(j,i)-1)))*T2(j,i);
        D3(j,i)=((M1(1,25))*sin((pi/5)*sector(j,i)-angle(j,i)))*T3(j,i);
        D4(j,i)=((-M1(1,16))*sin(angle(j,i)-(pi/5)*(sector(j,i)-1)))*T4(j,i);
        Dz(j,i)=0;
    end
end
end
for j=1:1:r
for i=1:fsw*time%Plane1
    P1 = 0;%HDF012347
    P2 = Qz(j,i)/2;
    P3 = Qz(j,i)/2+Q1(j,i);
    P4 = Qz(j,i)/2+Q1(j,i)+Q2(j,i);
    P5 = Qz(j,i)/2+Q1(j,i)+Q2(j,i)+Q3(j,i);
    P6 = Qz(j,i)/2+Q1(j,i)+Q2(j,i)+Q3(j,i)+Q4(j,i);
    P7 = Qz(j,i)+Q1(j,i)+Q2(j,i)+Q3(j,i)+Q4(j,i);
    HDFQ_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*Tz(j,i)/2;
    R1 = 0;
    R2 = Dz(j,i)/2;
    R3 = Dz(j,i)/2+D1(j,i);
    R4 = Dz(j,i)/2+D1(j,i)+D2(j,i);
    R5 = Dz(j,i)/2+D1(j,i)+D2(j,i)+D3(j,i);
    R6 = Dz(j,i)/2+D1(j,i)+D2(j,i)+D3(j,i)+D4(j,i);
    R7 = Dz(j,i)+D1(j,i)+D2(j,i)+D3(j,i)+D4(j,i);
    HDFD_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*Tz(j,i)/2;
    P1 = 0;%HDF012343
    P2 = Qz(j,i);
    P3 = Qz(j,i)+Q1(j,i);
    P4 = Qz(j,i)+Q1(j,i)+Q2(j,i);
    P5 = Qz(j,i)+Q1(j,i)+Q2(j,i)+Q3(j,i)*fraction3;
    P6 = Qz(j,i)+Q1(j,i)+Q2(j,i)+Q3(j,i)*fraction3+Q4(j,i);
    P7 = Qz(j,i)+Q1(j,i)+Q2(j,i)+Q3(j,i)+Q4(j,i);
    HDF0121Q_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Dz(j,i);
    R3 = Dz(j,i)+D1(j,i);
    R4 = Dz(j,i)+D1(j,i)+D2(j,i);
    R5 = Dz(j,i)+D1(j,i)+D2(j,i)+D3(j,i)*fraction3;
    R6 = Dz(j,i)+D1(j,i)+D2(j,i)+D3(j,i)*fraction3+D4(j,i);
    R7 = Dz(j,i)+D1(j,i)+D2(j,i)+D3(j,i)+D4(j,i);
    HDF0121D_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(j,i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Qz(j,i);
    P3 = Qz(j,i)+Q4(j,i);
    P4 = Qz(j,i)+Q4(j,i)+Q3(j,i);
    P5 = Qz(j,i)+Q4(j,i)+Q3(j,i)+Q2(j,i)*fraction3;
    P6 = Qz(j,i)+Q4(j,i)+Q3(j,i)+Q2(j,i)*fraction3+Q1(j,i);
    P7 = Qz(j,i)+Q4(j,i)+Q3(j,i)+Q2(j,i)+Q1(j,i);
    HDF7212Q_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Dz(j,i);
    R3 = Dz(j,i)+D4(j,i);
    R4 = Dz(j,i)+D4(j,i)+D3(j,i);
    R5 = Dz(j,i)+D4(j,i)+D3(j,i)+D2(j,i)*fraction3;
    R6 = Dz(j,i)+D4(j,i)+D3(j,i)+D2(j,i)*fraction3+D1(j,i);
    R7 = Dz(j,i)+D4(j,i)+D3(j,i)+D2(j,i)+D1(j,i);
    HDF7212D_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*(1-fraction3);
end
end
for j=1:1:r
for i=1:fsw*time%Plane2
    P1 = 0;%HDF012347
    P2 = Q9(j,i)/2;
    P3 = Q9(j,i)/2+Q5(j,i);
    P4 = Q9(j,i)/2+Q5(j,i)+Q6(j,i);
    P5 = Q9(j,i)/2+Q5(j,i)+Q6(j,i)+Q7(j,i);
    P6 = Q9(j,i)/2+Q5(j,i)+Q6(j,i)+Q7(j,i)+Q8(j,i);
    P7 = Q9(j,i)+Q5(j,i)+Q6(j,i)+Q7(j,i)+Q8(j,i);
    HDFQ_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*Tz(j,i)/2;
    R1 = 0;
    R2 = D9(j,i)/2;
    R3 = D9(j,i)/2+D5(j,i);
    R4 = D9(j,i)/2+D5(j,i)+D6(j,i);
    R5 = D9(j,i)/2+D5(j,i)+D6(j,i)+D7(j,i);
    R6 = D9(j,i)/2+D5(j,i)+D6(j,i)+D7(j,i)+D8(j,i);
    R7 = D9(j,i)+D5(j,i)+D6(j,i)+D7(j,i)+D8(j,i);
    HDFD_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*Tz(j,i)/2;
    P1 = 0;%HDF012343
    P2 = Q9(j,i);
    P3 = Q9(j,i)+Q5(j,i);
    P4 = Q9(j,i)+Q5(j,i)+Q6(j,i);
    P5 = Q9(j,i)+Q5(j,i)+Q6(j,i)+Q7(j,i)*fraction3;
    P6 = Q9(j,i)+Q5(j,i)+Q6(j,i)+Q7(j,i)*fraction3+Q8(j,i);
    P7 = Q9(j,i)+Q5(j,i)+Q6(j,i)+Q7(j,i)+Q8(j,i);
    HDF0121Q_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(j,i)*(1-fraction3);
    R1 = 0;
    R2 = D9(j,i);
    R3 = D9(j,i)+D5(j,i);
    R4 = D9(j,i)+D5(j,i)+D6(j,i);
    R5 = D9(j,i)+D5(j,i)+D6(j,i)+D7(j,i)*fraction3;
    R6 = D9(j,i)+D5(j,i)+D6(j,i)+D7(j,i)*fraction3+D8(j,i);
    R7 = D9(j,i)+D5(j,i)+D6(j,i)+D7(j,i)+D8(j,i);
    HDF0121D_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(j,i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Q9(j,i);
    P3 = Q9(j,i)+Q8(j,i);
    P4 = Q9(j,i)+Q8(j,i)+Q7(j,i);
    P5 = Q9(j,i)+Q8(j,i)+Q7(j,i)+Q6(j,i)*fraction3;
    P6 = Q9(j,i)+Q8(j,i)+Q7(j,i)+Q6(j,i)*fraction3+Q5(j,i);
    P7 = Q9(j,i)+Q8(j,i)+Q7(j,i)+Q6(j,i)+Q5(j,i);
    HDF7212Q_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*(1-fraction3);
    R1 = 0;
    R2 = D9(j,i);
    R3 = D9(j,i)+D8(j,i);
    R4 = D9(j,i)+D8(j,i)+D7(j,i);
    R5 = D9(j,i)+D8(j,i)+D7(j,i)+D6(j,i)*fraction3;
    R6 = D9(j,i)+D8(j,i)+D7(j,i)+D6(j,i)*fraction3+D5(j,i);
    R7 = D9(j,i)+D8(j,i)+D7(j,i)+D6(j,i)+D5(j,i);
    HDF7212D_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*(1-fraction3);
end
end
HDF = HDFQ_p1 + HDFD_p1 + HDFQ_p2 + HDFD_p2;
HDF0121 = HDF0121Q_p1 + HDF0121D_p1 + HDF0121Q_p2 + HDF0121D_p2;
HDF7212 = HDF7212Q_p1 + HDF7212D_p1 + HDF7212Q_p2 + HDF7212D_p2;

threeHDF=zeros(3*r,fsw*time);
minHDFcase=zeros(r,fsw*time);
minHDFvalue=zeros(r,fsw*time);
for j=1:1:r
    threeHDF(3*j-2,:)=HDF(j,:);
    threeHDF(3*j-1,:)=HDF0121(j,:);
    threeHDF(3*j,:)=HDF7212(j,:);
    minHDFvalue(j,:)=min(threeHDF(3*j-2:3*j,:));
end
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
sinHDFdistribute1x=zeros(r,fsw*time);
sinHDFdistribute1y=zeros(r,fsw*time);
sinHDFdistribute2x=zeros(r,fsw*time);
sinHDFdistribute2y=zeros(r,fsw*time);
sinHDFdistribute3x=zeros(r,fsw*time);
sinHDFdistribute3y=zeros(r,fsw*time);
figure;
for j=1:1:r
    plot(sindq_p1x(j,:),sindq_p1y(j,:),'b.');
    hold on;
end
 title('sinwave');
 axis([-0.6 0.6 -0.6 0.6],'square')
figure;%HDF蛾だガ
for j=1:1:r
for i=1:1:fsw*time
    switch minHDFcase(j,i)
        case 1
        sinHDFdistribute1x(j,i)=sindq_p1x(j,i);
        sinHDFdistribute1y(j,i)=sindq_p1y(j,i);
        case 2
        sinHDFdistribute2x(j,i)=sindq_p1x(j,i);
        sinHDFdistribute2y(j,i)=sindq_p1y(j,i);
        case 3
        sinHDFdistribute3x(j,i)=sindq_p1x(j,i);
        sinHDFdistribute3y(j,i)=sindq_p1y(j,i);
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
figure,plot(0:interval:amp,Fdist,'r');
hold on,plot(0:interval:amp,Fdist0121,'b');
hold on,plot(0:interval:amp,Fdist7212,'g');
title('Fdist');

figure;%HDF
plot(angle(r,:)*180/pi,HDF(r,:),'r');
hold on;
plot(angle(r,:)*180/pi,HDF0121(r,:),'b');
hold on;
plot(angle(r,:)*180/pi,HDF7212(r,:),'g');
title('HDF');

% figure;
% plot(angle(r,:)*180/pi,Q1(r,:),'r');
% hold on;
% plot(angle(r,:)*180/pi,Q2(r,:),'b');
% hold on;
% plot(angle(r,:)*180/pi,Q3(r,:),'g');
% hold on;
% plot(angle(r,:)*180/pi,Q4(r,:),'y');

% figure;%the same
% plot(angle(r,:)*180/pi,T1(r,:),'r');
% hold on;
% plot(angle(r,:)*180/pi,T2(r,:),'b');
% hold on;
% plot(angle(r,:)*180/pi,T3(r,:),'g');
% hold on;
% plot(angle(r,:)*180/pi,T4(r,:),'y');

% figure;
% plot(angle(r,:)*180/pi,sector(r,:),'b');

% figure;
% plot(angle(r,:),'b');
%}