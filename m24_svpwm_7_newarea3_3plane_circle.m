clear%5/17穝跋
close all
clc
start_time=clock;
time = 0.02;
fsw = 15000;
bit = 400;
amp = 0.53;%0.53
fraction3=0.5;
f = 50;
t = 0:1/fsw:time-1/fsw; 
interval = 0.01;
r = fix(amp/interval);
for j=1:1:r
    amp = interval*(j-1);
    s1(j,:) = amp*cos(2*pi*f*t);
    s2(j,:) = amp*cos(2*pi*f*t-2*pi/7); 
    s3(j,:) = amp*cos(2*pi*f*t-4*pi/7);
    s4(j,:) = amp*cos(2*pi*f*t-6*pi/7); 
    s5(j,:) = amp*cos(2*pi*f*t-8*pi/7); 
    s6(j,:) = amp*cos(2*pi*f*t-10*pi/7); 
    s7(j,:) = amp*cos(2*pi*f*t-12*pi/7); 
end
a = 2*pi/7;
mapping = 2/7*[...
    cos(0) cos(a) cos(2*a) cos(3*a) cos(4*a) cos(5*a) cos(6*a);...
    sin(0) sin(a) sin(2*a) sin(3*a) sin(4*a) sin(5*a) sin(6*a);...
    cos(0) cos(2*a) cos(4*a) cos(6*a) cos(8*a) cos(10*a) cos(12*a);...
    sin(0) sin(2*a) sin(4*a) sin(6*a) sin(8*a) sin(10*a) sin(12*a);...
    cos(0) cos(3*a) cos(6*a) cos(9*a) cos(12*a) cos(15*a) cos(18*a);...
    sin(0) sin(3*a) sin(6*a) sin(9*a) sin(12*a) sin(15*a) sin(18*a);...
    1 1 1 1 1 1 1];
for j=1:1:r
    sindq(7*j-6,:) = mapping(1,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq(7*j-5,:) = mapping(2,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq(7*j-4,:) = mapping(3,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq(7*j-3,:) = mapping(4,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq(7*j-2,:) = mapping(5,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq(7*j-1,:) = mapping(6,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq(7*j,:) = mapping(7,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:);s6(j,:);s7(j,:)];
    sindq_p1x(j,:) = sindq(7*j-6,:);
    sindq_p1y(j,:) = sindq(7*j-5,:);
    sindq_p2x(j,:) = sindq(7*j-4,:);
    sindq_p2y(j,:) = sindq(7*j-3,:);
    sindq_p3x(j,:) = sindq(7*j-2,:);
    sindq_p3y(j,:) = sindq(7*j-1,:);
end
angle=atan2(sindq_p1y,sindq_p1x);
for i=1:1:fsw*time
    if angle(2,i)<0
        angle(2,i)=angle(2,i)+2*pi;
    end
end
for j=3:1:r
    angle(j,:)=angle(2,:);
end
% angle=angle+1e-10;
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
for j=1:1:r
for i=1:fsw*time%耞sector
    if angle(j,i)>0 && angle(j,i)<=pi/7;%sector1
        sector(j,i)=1;
        sb1(i)=64;
        sb2(i)=96;
        sb3(i)=97;
        sb4(i)=113;
        sb5(i)=115;
        sb6(i)=123;
    end
    if angle(j,i)>pi/7 && angle(j,i)<=2*pi/7;%sector2
        sector(j,i)=2;
        sb1(i)=32;
        sb2(i)=96;
        sb3(i)=112;
        sb4(i)=113;
        sb5(i)=121;
        sb6(i)=123;
    end   
    if angle(j,i)>2*pi/7 && angle(j,i)<=3*pi/7;%sector3
        sector(j,i)=3;
        sb1(i)=32;
        sb2(i)=48;
        sb3(i)=112;
        sb4(i)=120;
        sb5(i)=121;
        sb6(i)=125;
    end
    if angle(j,i)>3*pi/7 && angle(j,i)<=4*pi/7;%sector4
        sector(j,i)=4;
        sb1(i)=16;
        sb2(i)=48;
        sb3(i)=56;
        sb4(i)=120;
        sb5(i)=124;
        sb6(i)=125;
    end
    if angle(j,i)>4*pi/7 && angle(j,i)<=5*pi/7;%sector5
        sector(j,i)=5;
        sb1(i)=16;
        sb2(i)=24;
        sb3(i)=56;
        sb4(i)=60;
        sb5(i)=124;
        sb6(i)=126;
    end 
    if angle(j,i)>5*pi/7 && angle(j,i)<=6*pi/7;%sector6
        sector(j,i)=6;
        sb1(i)=8;
        sb2(i)=24;
        sb3(i)=28;
        sb4(i)=60;
        sb5(i)=62;
        sb6(i)=126;
    end
    if angle(j,i)>6*pi/7 && angle(j,i)<=pi;%sector7
        sector(j,i)=7;
        sb1(i)=8;
        sb2(i)=12;
        sb3(i)=28;
        sb4(i)=30;
        sb5(i)=62;
        sb6(i)=63;
    end
    if angle(j,i)>pi && angle(j,i)<=8*pi/7;%sector8
        sector(j,i)=8;
        sb1(i)=4;
        sb2(i)=12;
        sb3(i)=14;
        sb4(i)=30;
        sb5(i)=31;
        sb6(i)=63;
    end
    if angle(j,i)>8*pi/7 && angle(j,i)<=9*pi/7;%sector9
        sector(j,i)=9;
        sb1(i)=4;
        sb2(i)=6;
        sb3(i)=14;
        sb4(i)=15;
        sb5(i)=31;
        sb6(i)=95;
    end
    if angle(j,i)>9*pi/7 && angle(j,i)<=10*pi/7;%sector10
        sector(j,i)=10;
        sb1(i)=2;
        sb2(i)=6;
        sb3(i)=7;
        sb4(i)=15;
        sb5(i)=79;
        sb6(i)=95;
    end
    if angle(j,i)>10*pi/7 && angle(j,i)<=11*pi/7;%sector11
        sector(j,i)=11;
        sb1(i)=2;
        sb2(i)=3;
        sb3(i)=7;
        sb4(i)=71;
        sb5(i)=79;
        sb6(i)=111;
    end
    if angle(j,i)>11*pi/7 && angle(j,i)<=12*pi/7;%sector12
        sector(j,i)=12;
        sb1(i)=1;
        sb2(i)=3;
        sb3(i)=67;
        sb4(i)=71;
        sb5(i)=103;
        sb6(i)=111;
    end
    if angle(j,i)>12*pi/7 && angle(j,i)<=13*pi/7;%sector13
        sector(j,i)=13;
        sb1(i)=1;
        sb2(i)=65;
        sb3(i)=67;
        sb4(i)=99;
        sb5(i)=103;
        sb6(i)=119;
    end
    if angle(j,i)>13*pi/7 && angle(j,i)<=2*pi;%sector14
        sector(j,i)=14;
        sb1(i)=64;
        sb2(i)=65;
        sb3(i)=97;
        sb4(i)=99;
        sb5(i)=115;
        sb6(i)=119;
    end
end
end
for j=1:1:r
for i=1:fsw*time%Plane2,3 Q1,Q2.Qz
     duty(:,i)=[M1(:,sb1(i)) M1(:,sb2(i)) M1(:,sb3(i)) M1(:,sb4(i)) M1(:,sb5(i)) M1(:,sb6(i));...
               M2(:,sb1(i)) M2(:,sb2(i)) M2(:,sb3(i)) M2(:,sb4(i)) M2(:,sb5(i)) M2(:,sb6(i));...
               M3(:,sb1(i)) M3(:,sb2(i)) M3(:,sb3(i)) M3(:,sb4(i)) M3(:,sb5(i)) M3(:,sb6(i))]\...
              [sindq_p1x(j,i);sindq_p1y(j,i);sindq_p2x(j,i);sindq_p2y(j,i);sindq_p3x(j,i);sindq_p3y(j,i)];
    Q1_p2(j,i)=M2(1,sb1(i))*duty(1,i);
    Q2_p2(j,i)=M2(1,sb2(i))*duty(2,i);
    Q3_p2(j,i)=M2(1,sb3(i))*duty(3,i);
    Q4_p2(j,i)=M2(1,sb4(i))*duty(4,i);
    Q5_p2(j,i)=M2(1,sb5(i))*duty(5,i);
    Q6_p2(j,i)=M2(1,sb6(i))*duty(6,i);
    Qz_p2(j,i)=0;
    D1_p2(j,i)=M2(2,sb1(i))*duty(1,i);
    D2_p2(j,i)=M2(2,sb2(i))*duty(2,i);
    D3_p2(j,i)=M2(2,sb3(i))*duty(3,i);
    D4_p2(j,i)=M2(2,sb4(i))*duty(4,i);
    D5_p2(j,i)=M2(2,sb5(i))*duty(5,i);
    D6_p2(j,i)=M2(2,sb6(i))*duty(6,i);
    Dz_p2(j,i)=0;
    Q1_p3(j,i)=M3(1,sb1(i))*duty(1,i);
    Q2_p3(j,i)=M3(1,sb2(i))*duty(2,i);
    Q3_p3(j,i)=M3(1,sb3(i))*duty(3,i);
    Q4_p3(j,i)=M3(1,sb4(i))*duty(4,i);
    Q5_p3(j,i)=M3(1,sb5(i))*duty(5,i);
    Q6_p3(j,i)=M3(1,sb6(i))*duty(6,i);
    Qz_p3(j,i)=0;
    D1_p3(j,i)=M3(2,sb1(i))*duty(1,i);
    D2_p3(j,i)=M3(2,sb2(i))*duty(2,i);
    D3_p3(j,i)=M3(2,sb3(i))*duty(3,i);
    D4_p3(j,i)=M3(2,sb4(i))*duty(4,i);
    D5_p3(j,i)=M3(2,sb5(i))*duty(5,i);
    D6_p3(j,i)=M3(2,sb6(i))*duty(6,i);
    Dz_p3(j,i)=0;
    a1(j,i)=duty(1,i);
    a2(j,i)=duty(2,i);
    a3(j,i)=duty(3,i);
    a4(j,i)=duty(4,i);
    a5(j,i)=duty(5,i);
    a6(j,i)=duty(6,i);
end
end
az=1-a1-a2-a3-a4-a5-a6;

T1=a1;%礚&Τ弘非
T2=a2;
T3=a3;
T4=a4;
T5=a5;
T6=a6;
Tz=az;
for j=1:1:r
    amp = interval*(j-1);
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(j,i),2)==1%1,sector碞琌计0,sector碞琌案计
        Q1_p1(j,i)=((M1(1,64))*cos(angle(j,i)-(pi/7)*(sector(j,i)-1))-amp)*T1(j,i);
        Q2_p1(j,i)=((M1(1,115))*cos((pi/7)*sector(j,i)-angle(j,i))-amp)*T2(j,i);
        Q3_p1(j,i)=((M1(1,97))*cos(angle(j,i)-(pi/7)*(sector(j,i)-1))-amp)*T3(j,i);
        Q4_p1(j,i)=((M1(1,97))*cos((pi/7)*sector(j,i)-angle(j,i))-amp)*T4(j,i);
        Q5_p1(j,i)=((M1(1,115))*cos(angle(j,i)-(pi/7)*(sector(j,i)-1))-amp)*T5(j,i);
        Q6_p1(j,i)=((M1(1,64))*cos((pi/7)*sector(j,i)-angle(j,i))-amp)*T6(j,i);
        Qz_p1(j,i)=-amp*Tz(j,i);
        D1_p1(j,i)=((M1(1,64))*sin(angle(j,i)-(pi/7)*(sector(j,i)-1)))*T1(j,i);
        D2_p1(j,i)=((-M1(1,115))*sin((pi/7)*sector(j,i)-angle(j,i)))*T2(j,i);
        D3_p1(j,i)=((M1(1,97))*sin(angle(j,i)-(pi/7)*(sector(j,i)-1)))*T3(j,i);
        D4_p1(j,i)=((-M1(1,97))*sin((pi/7)*sector(j,i)-angle(j,i)))*T4(j,i);
        D5_p1(j,i)=((M1(1,115))*sin(angle(j,i)-(pi/7)*(sector(j,i)-1)))*T5(j,i);
        D6_p1(j,i)=((-M1(1,64))*sin((pi/7)*sector(j,i)-angle(j,i)))*T6(j,i);
        Dz_p1(j,i)=0;
    else
        Q1_p1(j,i)=((M1(1,64))*cos((pi/7)*sector(j,i)-angle(j,i))-amp)*T1(j,i);
        Q2_p1(j,i)=((M1(1,115))*cos(angle(j,i)-(pi/7)*(sector(j,i)-1))-amp)*T2(j,i);
        Q3_p1(j,i)=((M1(1,97))*cos((pi/7)*sector(j,i)-angle(j,i))-amp)*T3(j,i);
        Q4_p1(j,i)=((M1(1,97))*cos(angle(j,i)-(pi/7)*(sector(j,i)-1))-amp)*T4(j,i);
        Q5_p1(j,i)=((M1(1,115))*cos((pi/7)*sector(j,i)-angle(j,i))-amp)*T5(j,i);
        Q6_p1(j,i)=((M1(1,64))*cos(angle(j,i)-(pi/7)*(sector(j,i)-1))-amp)*T6(j,i);
        Qz_p1(j,i)=-amp*Tz(j,i); 
        D1_p1(j,i)=((M1(1,64))*sin((pi/7)*sector(j,i)-angle(j,i)))*T1(j,i);
        D2_p1(j,i)=((-M1(1,115))*sin(angle(j,i)-(pi/7)*(sector(j,i)-1)))*T2(j,i);
        D3_p1(j,i)=((M1(1,97))*sin((pi/7)*sector(j,i)-angle(j,i)))*T3(j,i);
        D4_p1(j,i)=((-M1(1,97))*sin(angle(j,i)-(pi/7)*(sector(j,i)-1)))*T4(j,i);
        D5_p1(j,i)=((M1(1,115))*sin((pi/7)*sector(j,i)-angle(j,i)))*T5(j,i);
        D6_p1(j,i)=((-M1(1,64))*sin(angle(j,i)-(pi/7)*(sector(j,i)-1)))*T6(j,i);
        Dz_p1(j,i)=0;
    end
end
end
QQQ_p1=Qz_p1+Q1_p1+Q2_p1+Q3_p1+Q4_p1+Q5_p1+Q6_p1;
DDD_p2=Dz_p1+D1_p1+D2_p1+D3_p1+D4_p1+D5_p1+D6_p1;
for j=1:1:r
for i=1:1:fsw*time%Plane1
    P1=0;%01234567
    P2=Qz_p1(j,i)/2;
    P3=Qz_p1(j,i)/2+Q1_p1(j,i);
    P4=Qz_p1(j,i)/2+Q1_p1(j,i)+Q2_p1(j,i);
    P5=Qz_p1(j,i)/2+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i);
    P6=Qz_p1(j,i)/2+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i);
    P7=Qz_p1(j,i)/2+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i)+Q5_p1(j,i);
    P8=Qz_p1(j,i)/2+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i)+Q5_p1(j,i)+Q6_p1(j,i);
    P9=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i)+Q5_p1(j,i)+Q6_p1(j,i);
    HDFQ_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i)...
              +((P7)^2+(P7)*(P8)+(P8)^2)*T6(j,i)...
              +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(j,i)/2;
    R1=0;
    R2=Dz_p1(j,i)/2;
    R3=Dz_p1(j,i)/2+D1_p1(j,i);
    R4=Dz_p1(j,i)/2+D1_p1(j,i)+D2_p1(j,i);
    R5=Dz_p1(j,i)/2+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i);
    R6=Dz_p1(j,i)/2+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i);
    R7=Dz_p1(j,i)/2+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i)+D5_p1(j,i);
    R8=Dz_p1(j,i)/2+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i)+D5_p1(j,i)+D6_p1(j,i);
    R9=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i)+D5_p1(j,i)+D6_p1(j,i);
    HDFD_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i)...
              +((R7)^2+(R7)*(R8)+(R8)^2)*T6(j,i)...
              +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(j,i)/2;
    P1=0;%01234565
    P2=Qz_p1(j,i);
    P3=Qz_p1(j,i)+Q1_p1(j,i);
    P4=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i);
    P5=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i);
    P6=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i);
    P7=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i)+Q5_p1(j,i)*fraction3;
    P8=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i)+Q5_p1(j,i)*fraction3+Q6_p1(j,i);
    P9=Qz_p1(j,i)+Q1_p1(j,i)+Q2_p1(j,i)+Q3_p1(j,i)+Q4_p1(j,i)+Q5_p1(j,i)+Q6_p1(j,i);
    HDF0121Q_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(j,i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(j,i)*(1-fraction3);
    R1=0;
    R2=Dz_p1(j,i);
    R3=Dz_p1(j,i)+D1_p1(j,i);
    R4=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i);
    R5=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i);
    R6=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i);
    R7=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i)+D5_p1(j,i)*fraction3;
    R8=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i)+D5_p1(j,i)*fraction3+D6_p1(j,i);
    R9=Dz_p1(j,i)+D1_p1(j,i)+D2_p1(j,i)+D3_p1(j,i)+D4_p1(j,i)+D5_p1(j,i)+D6_p1(j,i);
    HDF0121D_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(j,i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(j,i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz_p1(j,i);
    P3=Qz_p1(j,i)+Q6_p1(j,i);
    P4=Qz_p1(j,i)+Q6_p1(j,i)+Q5_p1(j,i);
    P5=Qz_p1(j,i)+Q6_p1(j,i)+Q5_p1(j,i)+Q4_p1(j,i);
    P6=Qz_p1(j,i)+Q6_p1(j,i)+Q5_p1(j,i)+Q4_p1(j,i)+Q3_p1(j,i);
    P7=Qz_p1(j,i)+Q6_p1(j,i)+Q5_p1(j,i)+Q4_p1(j,i)+Q3_p1(j,i)+Q2_p1(j,i)*fraction3;
    P8=Qz_p1(j,i)+Q6_p1(j,i)+Q5_p1(j,i)+Q4_p1(j,i)+Q3_p1(j,i)+Q2_p1(j,i)*fraction3+Q1_p1(j,i);
    P9=Qz_p1(j,i)+Q6_p1(j,i)+Q5_p1(j,i)+Q4_p1(j,i)+Q3_p1(j,i)+Q2_p1(j,i)+Q1_p1(j,i);
    HDF7212Q_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(j,i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(j,i)...
                  +((P8)^2+(P8)*(P1)+(P1)^2)*T2(j,i)*(1-fraction3);
    R1=0;
    R2=Dz_p1(j,i);
    R3=Dz_p1(j,i)+D6_p1(j,i);
    R4=Dz_p1(j,i)+D6_p1(j,i)+D5_p1(j,i);
    R5=Dz_p1(j,i)+D6_p1(j,i)+D5_p1(j,i)+D4_p1(j,i);
    R6=Dz_p1(j,i)+D6_p1(j,i)+D5_p1(j,i)+D4_p1(j,i)+D3_p1(j,i);
    R7=Dz_p1(j,i)+D6_p1(j,i)+D5_p1(j,i)+D4_p1(j,i)+D3_p1(j,i)+D2_p1(j,i)*fraction3;
    R8=Dz_p1(j,i)+D6_p1(j,i)+D5_p1(j,i)+D4_p1(j,i)+D3_p1(j,i)+D2_p1(j,i)*fraction3+D1_p1(j,i);
    R9=Dz_p1(j,i)+D6_p1(j,i)+D5_p1(j,i)+D4_p1(j,i)+D3_p1(j,i)+D2_p1(j,i)+D1_p1(j,i);
    HDF7212D_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(j,i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(j,i)...
                  +((R8)^2+(R8)*(R1)+(R1)^2)*T2(j,i)*(1-fraction3);
end
end
for j=1:1:r
for i=1:1:fsw*time%Plane2
    P1=0;%01234567
    P2=Qz_p2(j,i)/2;
    P3=Qz_p2(j,i)/2+Q1_p2(j,i);
    P4=Qz_p2(j,i)/2+Q1_p2(j,i)+Q2_p2(j,i);
    P5=Qz_p2(j,i)/2+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i);
    P6=Qz_p2(j,i)/2+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i);
    P7=Qz_p2(j,i)/2+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i)+Q5_p2(j,i);
    P8=Qz_p2(j,i)/2+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i)+Q5_p2(j,i)+Q6_p2(j,i);
    P9=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i)+Q5_p2(j,i)+Q6_p2(j,i);
    HDFQ_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i)...
              +((P7)^2+(P7)*(P8)+(P8)^2)*T6(j,i)...
              +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(j,i)/2;
    R1=0;
    R2=Dz_p2(j,i)/2;
    R3=Dz_p2(j,i)/2+D1_p2(j,i);
    R4=Dz_p2(j,i)/2+D1_p2(j,i)+D2_p2(j,i);
    R5=Dz_p2(j,i)/2+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i);
    R6=Dz_p2(j,i)/2+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i);
    R7=Dz_p2(j,i)/2+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i)+D5_p2(j,i);
    R8=Dz_p2(j,i)/2+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i)+D5_p2(j,i)+D6_p2(j,i);
    R9=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i)+D5_p2(j,i)+D6_p2(j,i);
    HDFD_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i)...
              +((R7)^2+(R7)*(R8)+(R8)^2)*T6(j,i)...
              +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(j,i)/2;
    P1=0;%01234565
    P2=Qz_p2(j,i);
    P3=Qz_p2(j,i)+Q1_p2(j,i);
    P4=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i);
    P5=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i);
    P6=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i);
    P7=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i)+Q5_p2(j,i)*fraction3;
    P8=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i)+Q5_p2(j,i)*fraction3+Q6_p2(j,i);
    P9=Qz_p2(j,i)+Q1_p2(j,i)+Q2_p2(j,i)+Q3_p2(j,i)+Q4_p2(j,i)+Q5_p2(j,i)+Q6_p2(j,i);
    HDF0121Q_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(j,i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(j,i)*(1-fraction3);
    R1=0;
    R2=Dz_p2(j,i);
    R3=Dz_p2(j,i)+D1_p2(j,i);
    R4=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i);
    R5=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i);
    R6=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i);
    R7=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i)+D5_p2(j,i)*fraction3;
    R8=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i)+D5_p2(j,i)*fraction3+D6_p2(j,i);
    R9=Dz_p2(j,i)+D1_p2(j,i)+D2_p2(j,i)+D3_p2(j,i)+D4_p2(j,i)+D5_p2(j,i)+D6_p2(j,i);
    HDF0121D_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(j,i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(j,i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz_p2(j,i);
    P3=Qz_p2(j,i)+Q6_p2(j,i);
    P4=Qz_p2(j,i)+Q6_p2(j,i)+Q5_p2(j,i);
    P5=Qz_p2(j,i)+Q6_p2(j,i)+Q5_p2(j,i)+Q4_p2(j,i);
    P6=Qz_p2(j,i)+Q6_p2(j,i)+Q5_p2(j,i)+Q4_p2(j,i)+Q3_p2(j,i);
    P7=Qz_p2(j,i)+Q6_p2(j,i)+Q5_p2(j,i)+Q4_p2(j,i)+Q3_p2(j,i)+Q2_p2(j,i)*fraction3;
    P8=Qz_p2(j,i)+Q6_p2(j,i)+Q5_p2(j,i)+Q4_p2(j,i)+Q3_p2(j,i)+Q2_p2(j,i)*fraction3+Q1_p2(j,i);
    P9=Qz_p2(j,i)+Q6_p2(j,i)+Q5_p2(j,i)+Q4_p2(j,i)+Q3_p2(j,i)+Q2_p2(j,i)+Q1_p2(j,i);
    HDF7212Q_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(j,i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(j,i)...
                  +((P8)^2+(P8)*(P1)+(P1)^2)*T2(j,i)*(1-fraction3);
    R1=0;
    R2=Dz_p2(j,i);
    R3=Dz_p2(j,i)+D6_p2(j,i);
    R4=Dz_p2(j,i)+D6_p2(j,i)+D5_p2(j,i);
    R5=Dz_p2(j,i)+D6_p2(j,i)+D5_p2(j,i)+D4_p2(j,i);
    R6=Dz_p2(j,i)+D6_p2(j,i)+D5_p2(j,i)+D4_p2(j,i)+D3_p2(j,i);
    R7=Dz_p2(j,i)+D6_p2(j,i)+D5_p2(j,i)+D4_p2(j,i)+D3_p2(j,i)+D2_p2(j,i)*fraction3;
    R8=Dz_p2(j,i)+D6_p2(j,i)+D5_p2(j,i)+D4_p2(j,i)+D3_p2(j,i)+D2_p2(j,i)*fraction3+D1_p2(j,i);
    R9=Dz_p2(j,i)+D6_p2(j,i)+D5_p2(j,i)+D4_p2(j,i)+D3_p2(j,i)+D2_p2(j,i)+D1_p2(j,i);
    HDF7212D_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(j,i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(j,i)...
                  +((R8)^2+(R8)*(R1)+(R1)^2)*T2(j,i)*(1-fraction3);
end
end
for j=1:1:r
for i=1:1:fsw*time%Plane3
    P1=0;%01234567
    P2=Qz_p3(j,i)/2;
    P3=Qz_p3(j,i)/2+Q1_p3(j,i);
    P4=Qz_p3(j,i)/2+Q1_p3(j,i)+Q2_p3(j,i);
    P5=Qz_p3(j,i)/2+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i);
    P6=Qz_p3(j,i)/2+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i);
    P7=Qz_p3(j,i)/2+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i)+Q5_p3(j,i);
    P8=Qz_p3(j,i)/2+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i)+Q5_p3(j,i)+Q6_p3(j,i);
    P9=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i)+Q5_p3(j,i)+Q6_p3(j,i);
    HDFQ_p3(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i)...
              +((P7)^2+(P7)*(P8)+(P8)^2)*T6(j,i)...
              +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(j,i)/2;
    R1=0;
    R2=Dz_p3(j,i)/2;
    R3=Dz_p3(j,i)/2+D1_p3(j,i);
    R4=Dz_p3(j,i)/2+D1_p3(j,i)+D2_p3(j,i);
    R5=Dz_p3(j,i)/2+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i);
    R6=Dz_p3(j,i)/2+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i);
    R7=Dz_p3(j,i)/2+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i)+D5_p3(j,i);
    R8=Dz_p3(j,i)/2+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i)+D5_p3(j,i)+D6_p3(j,i);
    R9=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i)+D5_p3(j,i)+D6_p3(j,i);
    HDFD_p3(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i)...
              +((R7)^2+(R7)*(R8)+(R8)^2)*T6(j,i)...
              +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(j,i)/2;
    P1=0;%01234565
    P2=Qz_p3(j,i);
    P3=Qz_p3(j,i)+Q1_p3(j,i);
    P4=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i);
    P5=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i);
    P6=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i);
    P7=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i)+Q5_p3(j,i)*fraction3;
    P8=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i)+Q5_p3(j,i)*fraction3+Q6_p3(j,i);
    P9=Qz_p3(j,i)+Q1_p3(j,i)+Q2_p3(j,i)+Q3_p3(j,i)+Q4_p3(j,i)+Q5_p3(j,i)+Q6_p3(j,i);
    HDF0121Q_p3(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(j,i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(j,i)*(1-fraction3);
    R1=0;
    R2=Dz_p3(j,i);
    R3=Dz_p3(j,i)+D1_p3(j,i);
    R4=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i);
    R5=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i);
    R6=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i);
    R7=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i)+D5_p3(j,i)*fraction3;
    R8=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i)+D5_p3(j,i)*fraction3+D6_p3(j,i);
    R9=Dz_p3(j,i)+D1_p3(j,i)+D2_p3(j,i)+D3_p3(j,i)+D4_p3(j,i)+D5_p3(j,i)+D6_p3(j,i);
    HDF0121D_p3(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(j,i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(j,i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz_p3(j,i);
    P3=Qz_p3(j,i)+Q6_p3(j,i);
    P4=Qz_p3(j,i)+Q6_p3(j,i)+Q5_p3(j,i);
    P5=Qz_p3(j,i)+Q6_p3(j,i)+Q5_p3(j,i)+Q4_p3(j,i);
    P6=Qz_p3(j,i)+Q6_p3(j,i)+Q5_p3(j,i)+Q4_p3(j,i)+Q3_p3(j,i);
    P7=Qz_p3(j,i)+Q6_p3(j,i)+Q5_p3(j,i)+Q4_p3(j,i)+Q3_p3(j,i)+Q2_p3(j,i)*fraction3;
    P8=Qz_p3(j,i)+Q6_p3(j,i)+Q5_p3(j,i)+Q4_p3(j,i)+Q3_p3(j,i)+Q2_p3(j,i)*fraction3+Q1_p3(j,i);
    P9=Qz_p3(j,i)+Q6_p3(j,i)+Q5_p3(j,i)+Q4_p3(j,i)+Q3_p3(j,i)+Q2_p3(j,i)+Q1_p3(j,i);
    HDF7212Q_p3(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(j,i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(j,i)...
                  +((P8)^2+(P8)*(P1)+(P1)^2)*T2(j,i)*(1-fraction3);
    R1=0;
    R2=Dz_p3(j,i);
    R3=Dz_p3(j,i)+D6_p3(j,i);
    R4=Dz_p3(j,i)+D6_p3(j,i)+D5_p3(j,i);
    R5=Dz_p3(j,i)+D6_p3(j,i)+D5_p3(j,i)+D4_p3(j,i);
    R6=Dz_p3(j,i)+D6_p3(j,i)+D5_p3(j,i)+D4_p3(j,i)+D3_p3(j,i);
    R7=Dz_p3(j,i)+D6_p3(j,i)+D5_p3(j,i)+D4_p3(j,i)+D3_p3(j,i)+D2_p3(j,i)*fraction3;
    R8=Dz_p3(j,i)+D6_p3(j,i)+D5_p3(j,i)+D4_p3(j,i)+D3_p3(j,i)+D2_p3(j,i)*fraction3+D1_p3(j,i);
    R9=Dz_p3(j,i)+D6_p3(j,i)+D5_p3(j,i)+D4_p3(j,i)+D3_p3(j,i)+D2_p3(j,i)+D1_p3(j,i);
    HDF7212D_p3(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(j,i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(j,i)...
                  +((R8)^2+(R8)*(R1)+(R1)^2)*T2(j,i)*(1-fraction3);
end
end
HDF = HDFQ_p1 + HDFD_p1 + HDFQ_p2 + HDFD_p2 + HDFQ_p3 + HDFD_p3;
HDF0121 = HDF0121Q_p1 + HDF0121D_p1 + HDF0121Q_p2 + HDF0121D_p2 + HDF0121Q_p3 + HDF0121D_p3;
HDF7212 = HDF7212Q_p1 + HDF7212D_p1 + HDF7212Q_p2 + HDF7212D_p2 + HDF7212Q_p3 + HDF7212D_p3;

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

end_time=clock;
execution_time=end_time-start_time;

% figure;%allHDF
% plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
% plot(angle*180/pi,HDF0121,'b'),hold on;
% plot(angle*180/pi,HDF7212,'g'),hold on;
figure;%HDF
plot(angle(r,:)*180/pi,HDF(r,:),'r'),hold on;
plot(angle(r,:)*180/pi,HDF0121(r,:),'b'),hold on;
plot(angle(r,:)*180/pi,HDF7212(r,:),'g');
title('HDF');
%}