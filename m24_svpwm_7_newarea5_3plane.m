clear%11/19穝跋
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.51;%===============================================================
fraction3=0.1;
fraction5=0.9;%============================================================
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
        sb1(i)=64;
        sb2(i)=96;
        sb3(i)=97;
        sb4(i)=113;
        sb5(i)=115;
        sb6(i)=123;
        v1(:,i)=[1;0;0;0;0;0;0];
        v2(:,i)=[1;1;0;0;0;0;0];
        v3(:,i)=[1;1;0;0;0;0;1];
        v4(:,i)=[1;1;1;0;0;0;1];
        v5(:,i)=[1;1;1;0;0;1;1];
        v6(:,i)=[1;1;1;1;0;1;1];
    end
    if angle(i)>pi/7 && angle(i)<=2*pi/7;%sector2
        sector(i)=2;
        sb1(i)=32;
        sb2(i)=96;
        sb3(i)=112;
        sb4(i)=113;
        sb5(i)=121;
        sb6(i)=123;
        v1(:,i)=[0;1;0;0;0;0;0];
        v2(:,i)=[1;1;0;0;0;0;0];
        v3(:,i)=[1;1;1;0;0;0;0];
        v4(:,i)=[1;1;1;0;0;0;1];
        v5(:,i)=[1;1;1;1;0;0;1];
        v6(:,i)=[1;1;1;1;0;1;1];
    end   
    if angle(i)>2*pi/7 && angle(i)<=3*pi/7;%sector3
        sector(i)=3;
        sb1(i)=32;
        sb2(i)=48;
        sb3(i)=112;
        sb4(i)=120;
        sb5(i)=121;
        sb6(i)=125;
        v1(:,i)=[0;1;0;0;0;0;0];
        v2(:,i)=[0;1;1;0;0;0;0];
        v3(:,i)=[1;1;1;0;0;0;0];
        v4(:,i)=[1;1;1;1;0;0;0];
        v5(:,i)=[1;1;1;1;0;0;1];
        v6(:,i)=[1;1;1;1;1;0;1];
    end
    if angle(i)>3*pi/7 && angle(i)<=4*pi/7;%sector4
        sector(i)=4;
        sb1(i)=16;
        sb2(i)=48;
        sb3(i)=56;
        sb4(i)=120;
        sb5(i)=124;
        sb6(i)=125;
        v1(:,i)=[0;0;1;0;0;0;0];
        v2(:,i)=[0;1;1;0;0;0;0];
        v3(:,i)=[0;1;1;1;0;0;0];
        v4(:,i)=[1;1;1;1;0;0;0];
        v5(:,i)=[1;1;1;1;1;0;0];
        v6(:,i)=[1;1;1;1;1;0;1];
    end
    if angle(i)>4*pi/7 && angle(i)<=5*pi/7;%sector5
        sector(i)=5;
        sb1(i)=16;
        sb2(i)=24;
        sb3(i)=56;
        sb4(i)=60;
        sb5(i)=124;
        sb6(i)=126;
        v1(:,i)=[0;0;1;0;0;0;0];
        v2(:,i)=[0;0;1;1;0;0;0];
        v3(:,i)=[0;1;1;1;0;0;0];
        v4(:,i)=[0;1;1;1;1;0;0];
        v5(:,i)=[1;1;1;1;1;0;0];
        v6(:,i)=[1;1;1;1;1;1;0];
    end 
    if angle(i)>5*pi/7 && angle(i)<=6*pi/7;%sector6
        sector(i)=6;
        sb1(i)=8;
        sb2(i)=24;
        sb3(i)=28;
        sb4(i)=60;
        sb5(i)=62;
        sb6(i)=126;
        v1(:,i)=[0;0;0;1;0;0;0];
        v2(:,i)=[0;0;1;1;0;0;0];
        v3(:,i)=[0;0;1;1;1;0;0];
        v4(:,i)=[0;1;1;1;1;0;0];
        v5(:,i)=[0;1;1;1;1;1;0];
        v6(:,i)=[1;1;1;1;1;1;0];
    end
    if angle(i)>6*pi/7 && angle(i)<=pi;%sector7
        sector(i)=7;
        sb1(i)=8;
        sb2(i)=12;
        sb3(i)=28;
        sb4(i)=30;
        sb5(i)=62;
        sb6(i)=63;
        v1(:,i)=[0;0;0;1;0;0;0];
        v2(:,i)=[0;0;0;1;1;0;0];
        v3(:,i)=[0;0;1;1;1;0;0];
        v4(:,i)=[0;0;1;1;1;1;0];
        v5(:,i)=[0;1;1;1;1;1;0];
        v6(:,i)=[0;1;1;1;1;1;1];
    end
    if angle(i)>pi && angle(i)<=8*pi/7;%sector8
        sector(i)=8;
        sb1(i)=4;
        sb2(i)=12;
        sb3(i)=14;
        sb4(i)=30;
        sb5(i)=31;
        sb6(i)=63;
        v1(:,i)=[0;0;0;0;1;0;0];
        v2(:,i)=[0;0;0;1;1;0;0];
        v3(:,i)=[0;0;0;1;1;1;0];
        v4(:,i)=[0;0;1;1;1;1;0];
        v5(:,i)=[0;0;1;1;1;1;1];
        v6(:,i)=[0;1;1;1;1;1;1];
    end
    if angle(i)>8*pi/7 && angle(i)<=9*pi/7;%sector9
        sector(i)=9;
        sb1(i)=4;
        sb2(i)=6;
        sb3(i)=14;
        sb4(i)=15;
        sb5(i)=31;
        sb6(i)=95;
        v1(:,i)=[0;0;0;0;1;0;0];
        v2(:,i)=[0;0;0;0;1;1;0];
        v3(:,i)=[0;0;0;1;1;1;0];
        v4(:,i)=[0;0;0;1;1;1;1];
        v5(:,i)=[0;0;1;1;1;1;1];
        v6(:,i)=[1;0;1;1;1;1;1];
    end
    if angle(i)>9*pi/7 && angle(i)<=10*pi/7;%sector10
        sector(i)=10;
        sb1(i)=2;
        sb2(i)=6;
        sb3(i)=7;
        sb4(i)=15;
        sb5(i)=79;
        sb6(i)=95;
        v1(:,i)=[0;0;0;0;0;1;0];
        v2(:,i)=[0;0;0;0;1;1;0];
        v3(:,i)=[0;0;0;0;1;1;1];
        v4(:,i)=[0;0;0;1;1;1;1];
        v5(:,i)=[1;0;0;1;1;1;1];
        v6(:,i)=[1;0;1;1;1;1;1];
    end
    if angle(i)>10*pi/7 && angle(i)<=11*pi/7;%sector11
        sector(i)=11;
        sb1(i)=2;
        sb2(i)=3;
        sb3(i)=7;
        sb4(i)=71;
        sb5(i)=79;
        sb6(i)=111;
        v1(:,i)=[0;0;0;0;0;1;0];
        v2(:,i)=[0;0;0;0;0;1;1];
        v3(:,i)=[0;0;0;0;1;1;1];
        v4(:,i)=[1;0;0;0;1;1;1];
        v5(:,i)=[1;0;0;1;1;1;1];
        v6(:,i)=[1;1;0;1;1;1;1];
    end
    if angle(i)>11*pi/7 && angle(i)<=12*pi/7;%sector12
        sector(i)=12;
        sb1(i)=1;
        sb2(i)=3;
        sb3(i)=67;
        sb4(i)=71;
        sb5(i)=103;
        sb6(i)=111;
        v1(:,i)=[0;0;0;0;0;0;1];
        v2(:,i)=[0;0;0;0;0;1;1];
        v3(:,i)=[1;0;0;0;0;1;1];
        v4(:,i)=[1;0;0;0;1;1;1];
        v5(:,i)=[1;1;0;0;1;1;1];
        v6(:,i)=[1;1;0;1;1;1;1];
    end
    if angle(i)>12*pi/7 && angle(i)<=13*pi/7;%sector13
        sector(i)=13;
        sb1(i)=1;
        sb2(i)=65;
        sb3(i)=67;
        sb4(i)=99;
        sb5(i)=103;
        sb6(i)=119;
        v1(:,i)=[0;0;0;0;0;0;1];
        v2(:,i)=[1;0;0;0;0;0;1];
        v3(:,i)=[1;0;0;0;0;1;1];
        v4(:,i)=[1;1;0;0;0;1;1];
        v5(:,i)=[1;1;0;0;1;1;1];
        v6(:,i)=[1;1;1;0;1;1;1];
    end
    if angle(i)>13*pi/7 && angle(i)<=2*pi;%sector14
        sector(i)=14;
        sb1(i)=64;
        sb2(i)=65;
        sb3(i)=97;
        sb4(i)=99;
        sb5(i)=115;
        sb6(i)=119;
        v1(:,i)=[1;0;0;0;0;0;0];
        v2(:,i)=[1;0;0;0;0;0;1];
        v3(:,i)=[1;1;0;0;0;0;1];
        v4(:,i)=[1;1;0;0;0;1;1];
        v5(:,i)=[1;1;1;0;0;1;1];
        v6(:,i)=[1;1;1;0;1;1;1];
    end
end
for i=1:fsw*time%Plane2,3 Q1,Q2.Qz
    duty(:,i)=[M1(:,sb1(i)) M1(:,sb2(i)) M1(:,sb3(i)) M1(:,sb4(i)) M1(:,sb5(i)) M1(:,sb6(i));...
               M2(:,sb1(i)) M2(:,sb2(i)) M2(:,sb3(i)) M2(:,sb4(i)) M2(:,sb5(i)) M2(:,sb6(i));...
               M3(:,sb1(i)) M3(:,sb2(i)) M3(:,sb3(i)) M3(:,sb4(i)) M3(:,sb5(i)) M3(:,sb6(i))]\...
              [sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i);sindq3x(i);sindq3y(i)];
    Q1_p2(i)=M2(1,sb1(i))*duty(1,i);
    Q2_p2(i)=M2(1,sb2(i))*duty(2,i);
    Q3_p2(i)=M2(1,sb3(i))*duty(3,i);
    Q4_p2(i)=M2(1,sb4(i))*duty(4,i);
    Q5_p2(i)=M2(1,sb5(i))*duty(5,i);
    Q6_p2(i)=M2(1,sb6(i))*duty(6,i);
    Qz_p2(i)=0;
    D1_p2(i)=M2(2,sb1(i))*duty(1,i);
    D2_p2(i)=M2(2,sb2(i))*duty(2,i);
    D3_p2(i)=M2(2,sb3(i))*duty(3,i);
    D4_p2(i)=M2(2,sb4(i))*duty(4,i);
    D5_p2(i)=M2(2,sb5(i))*duty(5,i);
    D6_p2(i)=M2(2,sb6(i))*duty(6,i);
    Dz_p2(i)=0;
    Q1_p3(i)=M3(1,sb1(i))*duty(1,i);
    Q2_p3(i)=M3(1,sb2(i))*duty(2,i);
    Q3_p3(i)=M3(1,sb3(i))*duty(3,i);
    Q4_p3(i)=M3(1,sb4(i))*duty(4,i);
    Q5_p3(i)=M3(1,sb5(i))*duty(5,i);
    Q6_p3(i)=M3(1,sb6(i))*duty(6,i);
    Qz_p3(i)=0;
    D1_p3(i)=M3(2,sb1(i))*duty(1,i);
    D2_p3(i)=M3(2,sb2(i))*duty(2,i);
    D3_p3(i)=M3(2,sb3(i))*duty(3,i);
    D4_p3(i)=M3(2,sb4(i))*duty(4,i);
    D5_p3(i)=M3(2,sb5(i))*duty(5,i);
    D6_p3(i)=M3(2,sb6(i))*duty(6,i);
    Dz_p3(i)=0;
end
a1=duty(1,:);
a2=duty(2,:);
a3=duty(3,:);
a4=duty(4,:);
a5=duty(5,:);
a6=duty(6,:);
az=1-a1-a2-a3-a4-a5-a6;

T1=a1;%礚&Τ弘非
T2=a2;
T3=a3;
T4=a4;
T5=a5;
T6=a6;
Tz=az;
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(i),2)==1%1,sector碞琌计0,sector碞琌案计
        Q1_p1(i)=((M1(1,64))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T1(i);
        Q2_p1(i)=((M1(1,115))*cos((pi/7)*sector(i)-angle(i))-amp)*T2(i);
        Q3_p1(i)=((M1(1,97))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T3(i);
        Q4_p1(i)=((M1(1,97))*cos((pi/7)*sector(i)-angle(i))-amp)*T4(i);
        Q5_p1(i)=((M1(1,115))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T5(i);
        Q6_p1(i)=((M1(1,64))*cos((pi/7)*sector(i)-angle(i))-amp)*T6(i);
        Qz_p1(i)=-amp*Tz(i);
        D1_p1(i)=((M1(1,64))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T1(i);
        D2_p1(i)=((-M1(1,115))*sin((pi/7)*sector(i)-angle(i)))*T2(i);
        D3_p1(i)=((M1(1,97))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T3(i);
        D4_p1(i)=((-M1(1,97))*sin((pi/7)*sector(i)-angle(i)))*T4(i);
        D5_p1(i)=((M1(1,115))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T5(i);
        D6_p1(i)=((-M1(1,64))*sin((pi/7)*sector(i)-angle(i)))*T6(i);
        Dz_p1(i)=0;
    else
        Q1_p1(i)=((M1(1,64))*cos((pi/7)*sector(i)-angle(i))-amp)*T1(i);
        Q2_p1(i)=((M1(1,115))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T2(i);
        Q3_p1(i)=((M1(1,97))*cos((pi/7)*sector(i)-angle(i))-amp)*T3(i);
        Q4_p1(i)=((M1(1,97))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T4(i);
        Q5_p1(i)=((M1(1,115))*cos((pi/7)*sector(i)-angle(i))-amp)*T5(i);
        Q6_p1(i)=((M1(1,64))*cos(angle(i)-(pi/7)*(sector(i)-1))-amp)*T6(i);
        Qz_p1(i)=-amp*Tz(i); 
        D1_p1(i)=((M1(1,64))*sin((pi/7)*sector(i)-angle(i)))*T1(i);
        D2_p1(i)=((-M1(1,115))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T2(i);
        D3_p1(i)=((M1(1,97))*sin((pi/7)*sector(i)-angle(i)))*T3(i);
        D4_p1(i)=((-M1(1,97))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T4(i);
        D5_p1(i)=((M1(1,115))*sin((pi/7)*sector(i)-angle(i)))*T5(i);
        D6_p1(i)=((-M1(1,64))*sin(angle(i)-(pi/7)*(sector(i)-1)))*T6(i);
        Dz_p1(i)=0;
    end
end
QQQ_p1=Qz_p1+Q1_p1+Q2_p1+Q3_p1+Q4_p1+Q5_p1+Q6_p1;
DDD_p2=Dz_p1+D1_p1+D2_p1+D3_p1+D4_p1+D5_p1+D6_p1;
for i=1:1:fsw*time%Plane1
    P1=0;%01234567
    P2=Qz_p1(i)/2;
    P3=Qz_p1(i)/2+Q1_p1(i);
    P4=Qz_p1(i)/2+Q1_p1(i)+Q2_p1(i);
    P5=Qz_p1(i)/2+Q1_p1(i)+Q2_p1(i)+Q3_p1(i);
    P6=Qz_p1(i)/2+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i);
    P7=Qz_p1(i)/2+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i);
    P8=Qz_p1(i)/2+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)+Q6_p1(i);
    P9=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)+Q6_p1(i);
    HDFQ_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)...
              +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
              +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Dz_p1(i)/2;
    R3=Dz_p1(i)/2+D1_p1(i);
    R4=Dz_p1(i)/2+D1_p1(i)+D2_p1(i);
    R5=Dz_p1(i)/2+D1_p1(i)+D2_p1(i)+D3_p1(i);
    R6=Dz_p1(i)/2+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i);
    R7=Dz_p1(i)/2+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i);
    R8=Dz_p1(i)/2+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)+D6_p1(i);
    R9=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)+D6_p1(i);
    HDFD_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)...
              +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
              +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Qz_p1(i);
    P3=Qz_p1(i)+Q1_p1(i);
    P4=Qz_p1(i)+Q1_p1(i)+Q2_p1(i);
    P5=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i);
    P6=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i);
    P7=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)*fraction3;
    P8=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)*fraction3+Q6_p1(i);
    P9=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)+Q6_p1(i);
    HDF0121Q_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Dz_p1(i);
    R3=Dz_p1(i)+D1_p1(i);
    R4=Dz_p1(i)+D1_p1(i)+D2_p1(i);
    R5=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i);
    R6=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i);
    R7=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)*fraction3;
    R8=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)*fraction3+D6_p1(i);
    R9=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)+D6_p1(i);
    HDF0121D_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz_p1(i);
    P3=Qz_p1(i)+Q6_p1(i);
    P4=Qz_p1(i)+Q6_p1(i)+Q5_p1(i);
    P5=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i);
    P6=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i);
    P7=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i)+Q2_p1(i)*fraction3;
    P8=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i)+Q2_p1(i)*fraction3+Q1_p1(i);
    P9=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i)+Q2_p1(i)+Q1_p1(i);
    HDF7212Q_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                  +((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Dz_p1(i);
    R3=Dz_p1(i)+D6_p1(i);
    R4=Dz_p1(i)+D6_p1(i)+D5_p1(i);
    R5=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i);
    R6=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i);
    R7=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i)+D2_p1(i)*fraction3;
    R8=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i)+D2_p1(i)*fraction3+D1_p1(i);
    R9=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i)+D2_p1(i)+D1_p1(i);
    HDF7212D_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                  +((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction3);
    P1=0;%new01234565
    P2=Qz_p1(i);
    P3=Qz_p1(i)+Q1_p1(i);
    P4=Qz_p1(i)+Q1_p1(i)+Q2_p1(i);
    P5=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i);
    P6=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i);
    P7=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)*fraction5;
    P8=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)*fraction5+Q6_p1(i);
    P9=Qz_p1(i)+Q1_p1(i)+Q2_p1(i)+Q3_p1(i)+Q4_p1(i)+Q5_p1(i)+Q6_p1(i);
    HDFnew0121Q_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                     +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                     +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                     +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                     +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                     +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction5...
                     +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                     +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction5);
    R1=0;
    R2=Dz_p1(i);
    R3=Dz_p1(i)+D1_p1(i);
    R4=Dz_p1(i)+D1_p1(i)+D2_p1(i);
    R5=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i);
    R6=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i);
    R7=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)*fraction5;
    R8=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)*fraction5+D6_p1(i);
    R9=Dz_p1(i)+D1_p1(i)+D2_p1(i)+D3_p1(i)+D4_p1(i)+D5_p1(i)+D6_p1(i);
    HDFnew0121D_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                     +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                     +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                     +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                     +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                     +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction5...
                     +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                     +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction5);
    P1=0;%new76543212
    P2=Qz_p1(i);
    P3=Qz_p1(i)+Q6_p1(i);
    P4=Qz_p1(i)+Q6_p1(i)+Q5_p1(i);
    P5=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i);
    P6=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i);
    P7=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i)+Q2_p1(i)*fraction5;
    P8=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i)+Q2_p1(i)*fraction5+Q1_p1(i);
    P9=Qz_p1(i)+Q6_p1(i)+Q5_p1(i)+Q4_p1(i)+Q3_p1(i)+Q2_p1(i)+Q1_p1(i);
    HDFnew7212Q_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                     +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                     +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                     +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                     +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                     +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction5...
                     +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                     +((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction5);
    R1=0;
    R2=Dz_p1(i);
    R3=Dz_p1(i)+D6_p1(i);
    R4=Dz_p1(i)+D6_p1(i)+D5_p1(i);
    R5=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i);
    R6=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i);
    R7=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i)+D2_p1(i)*fraction5;
    R8=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i)+D2_p1(i)*fraction5+D1_p1(i);
    R9=Dz_p1(i)+D6_p1(i)+D5_p1(i)+D4_p1(i)+D3_p1(i)+D2_p1(i)+D1_p1(i);
    HDFnew7212D_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                     +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                     +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                     +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                     +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                     +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction5...
                     +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                     +((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction5);
end
for i=1:1:fsw*time%Plane2
    P1=0;%01234567
    P2=Qz_p2(i)/2;
    P3=Qz_p2(i)/2+Q1_p2(i);
    P4=Qz_p2(i)/2+Q1_p2(i)+Q2_p2(i);
    P5=Qz_p2(i)/2+Q1_p2(i)+Q2_p2(i)+Q3_p2(i);
    P6=Qz_p2(i)/2+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i);
    P7=Qz_p2(i)/2+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i);
    P8=Qz_p2(i)/2+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)+Q6_p2(i);
    P9=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)+Q6_p2(i);
    HDFQ_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)...
              +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
              +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Dz_p2(i)/2;
    R3=Dz_p2(i)/2+D1_p2(i);
    R4=Dz_p2(i)/2+D1_p2(i)+D2_p2(i);
    R5=Dz_p2(i)/2+D1_p2(i)+D2_p2(i)+D3_p2(i);
    R6=Dz_p2(i)/2+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i);
    R7=Dz_p2(i)/2+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i);
    R8=Dz_p2(i)/2+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)+D6_p2(i);
    R9=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)+D6_p2(i);
    HDFD_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)...
              +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
              +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Qz_p2(i);
    P3=Qz_p2(i)+Q1_p2(i);
    P4=Qz_p2(i)+Q1_p2(i)+Q2_p2(i);
    P5=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i);
    P6=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i);
    P7=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)*fraction3;
    P8=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)*fraction3+Q6_p2(i);
    P9=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)+Q6_p2(i);
    HDF0121Q_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Dz_p2(i);
    R3=Dz_p2(i)+D1_p2(i);
    R4=Dz_p2(i)+D1_p2(i)+D2_p2(i);
    R5=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i);
    R6=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i);
    R7=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)*fraction3;
    R8=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)*fraction3+D6_p2(i);
    R9=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)+D6_p2(i);
    HDF0121D_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz_p2(i);
    P3=Qz_p2(i)+Q6_p2(i);
    P4=Qz_p2(i)+Q6_p2(i)+Q5_p2(i);
    P5=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i);
    P6=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i);
    P7=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i)+Q2_p2(i)*fraction3;
    P8=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i)+Q2_p2(i)*fraction3+Q1_p2(i);
    P9=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i)+Q2_p2(i)+Q1_p2(i);
    HDF7212Q_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                  +((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Dz_p2(i);
    R3=Dz_p2(i)+D6_p2(i);
    R4=Dz_p2(i)+D6_p2(i)+D5_p2(i);
    R5=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i);
    R6=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i);
    R7=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i)+D2_p2(i)*fraction3;
    R8=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i)+D2_p2(i)*fraction3+D1_p2(i);
    R9=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i)+D2_p2(i)+D1_p2(i);
    HDF7212D_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                  +((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction3);
    P1=0;%new01234565
    P2=Qz_p2(i);
    P3=Qz_p2(i)+Q1_p2(i);
    P4=Qz_p2(i)+Q1_p2(i)+Q2_p2(i);
    P5=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i);
    P6=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i);
    P7=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)*fraction5;
    P8=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)*fraction5+Q6_p2(i);
    P9=Qz_p2(i)+Q1_p2(i)+Q2_p2(i)+Q3_p2(i)+Q4_p2(i)+Q5_p2(i)+Q6_p2(i);
    HDFnew0121Q_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                     +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                     +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                     +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                     +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                     +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction5...
                     +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                     +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction5);
    R1=0;
    R2=Dz_p2(i);
    R3=Dz_p2(i)+D1_p2(i);
    R4=Dz_p2(i)+D1_p2(i)+D2_p2(i);
    R5=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i);
    R6=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i);
    R7=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)*fraction5;
    R8=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)*fraction5+D6_p2(i);
    R9=Dz_p2(i)+D1_p2(i)+D2_p2(i)+D3_p2(i)+D4_p2(i)+D5_p2(i)+D6_p2(i);
    HDFnew0121D_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                     +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                     +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                     +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                     +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                     +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction5...
                     +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                     +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction5);
    P1=0;%new76543212
    P2=Qz_p2(i);
    P3=Qz_p2(i)+Q6_p2(i);
    P4=Qz_p2(i)+Q6_p2(i)+Q5_p2(i);
    P5=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i);
    P6=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i);
    P7=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i)+Q2_p2(i)*fraction5;
    P8=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i)+Q2_p2(i)*fraction5+Q1_p2(i);
    P9=Qz_p2(i)+Q6_p2(i)+Q5_p2(i)+Q4_p2(i)+Q3_p2(i)+Q2_p2(i)+Q1_p2(i);
    HDFnew7212Q_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                     +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                     +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                     +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                     +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                     +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction5...
                     +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                     +((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction5);
    R1=0;
    R2=Dz_p2(i);
    R3=Dz_p2(i)+D6_p2(i);
    R4=Dz_p2(i)+D6_p2(i)+D5_p2(i);
    R5=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i);
    R6=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i);
    R7=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i)+D2_p2(i)*fraction5;
    R8=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i)+D2_p2(i)*fraction5+D1_p2(i);
    R9=Dz_p2(i)+D6_p2(i)+D5_p2(i)+D4_p2(i)+D3_p2(i)+D2_p2(i)+D1_p2(i);
    HDFnew7212D_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                     +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                     +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                     +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                     +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                     +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction5...
                     +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                     +((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction5);
end
for i=1:1:fsw*time%Plane3
    P1=0;%01234567
    P2=Qz_p3(i)/2;
    P3=Qz_p3(i)/2+Q1_p3(i);
    P4=Qz_p3(i)/2+Q1_p3(i)+Q2_p3(i);
    P5=Qz_p3(i)/2+Q1_p3(i)+Q2_p3(i)+Q3_p3(i);
    P6=Qz_p3(i)/2+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i);
    P7=Qz_p3(i)/2+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i);
    P8=Qz_p3(i)/2+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)+Q6_p3(i);
    P9=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)+Q6_p3(i);
    HDFQ_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)...
              +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
              +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Dz_p3(i)/2;
    R3=Dz_p3(i)/2+D1_p3(i);
    R4=Dz_p3(i)/2+D1_p3(i)+D2_p3(i);
    R5=Dz_p3(i)/2+D1_p3(i)+D2_p3(i)+D3_p3(i);
    R6=Dz_p3(i)/2+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i);
    R7=Dz_p3(i)/2+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i);
    R8=Dz_p3(i)/2+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)+D6_p3(i);
    R9=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)+D6_p3(i);
    HDFD_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)...
              +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
              +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Qz_p3(i);
    P3=Qz_p3(i)+Q1_p3(i);
    P4=Qz_p3(i)+Q1_p3(i)+Q2_p3(i);
    P5=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i);
    P6=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i);
    P7=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)*fraction3;
    P8=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)*fraction3+Q6_p3(i);
    P9=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)+Q6_p3(i);
    HDF0121Q_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Dz_p3(i);
    R3=Dz_p3(i)+D1_p3(i);
    R4=Dz_p3(i)+D1_p3(i)+D2_p3(i);
    R5=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i);
    R6=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i);
    R7=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)*fraction3;
    R8=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)*fraction3+D6_p3(i);
    R9=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)+D6_p3(i);
    HDF0121D_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Qz_p3(i);
    P3=Qz_p3(i)+Q6_p3(i);
    P4=Qz_p3(i)+Q6_p3(i)+Q5_p3(i);
    P5=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i);
    P6=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i);
    P7=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i)+Q2_p3(i)*fraction3;
    P8=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i)+Q2_p3(i)*fraction3+Q1_p3(i);
    P9=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i)+Q2_p3(i)+Q1_p3(i);
    HDF7212Q_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                  +((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Dz_p3(i);
    R3=Dz_p3(i)+D6_p3(i);
    R4=Dz_p3(i)+D6_p3(i)+D5_p3(i);
    R5=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i);
    R6=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i);
    R7=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i)+D2_p3(i)*fraction3;
    R8=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i)+D2_p3(i)*fraction3+D1_p3(i);
    R9=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i)+D2_p3(i)+D1_p3(i);
    HDF7212D_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                  +((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction3);
    P1=0;%new01234565
    P2=Qz_p3(i);
    P3=Qz_p3(i)+Q1_p3(i);
    P4=Qz_p3(i)+Q1_p3(i)+Q2_p3(i);
    P5=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i);
    P6=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i);
    P7=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)*fraction5;
    P8=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)*fraction5+Q6_p3(i);
    P9=Qz_p3(i)+Q1_p3(i)+Q2_p3(i)+Q3_p3(i)+Q4_p3(i)+Q5_p3(i)+Q6_p3(i);
    HDFnew0121Q_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                     +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                     +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                     +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                     +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                     +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction5...
                     +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                     +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction5);
    R1=0;
    R2=Dz_p3(i);
    R3=Dz_p3(i)+D1_p3(i);
    R4=Dz_p3(i)+D1_p3(i)+D2_p3(i);
    R5=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i);
    R6=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i);
    R7=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)*fraction5;
    R8=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)*fraction5+D6_p3(i);
    R9=Dz_p3(i)+D1_p3(i)+D2_p3(i)+D3_p3(i)+D4_p3(i)+D5_p3(i)+D6_p3(i);
    HDFnew0121D_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                     +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                     +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                     +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                     +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                     +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction5...
                     +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                     +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction5);
    P1=0;%new76543212
    P2=Qz_p3(i);
    P3=Qz_p3(i)+Q6_p3(i);
    P4=Qz_p3(i)+Q6_p3(i)+Q5_p3(i);
    P5=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i);
    P6=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i);
    P7=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i)+Q2_p3(i)*fraction5;
    P8=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i)+Q2_p3(i)*fraction5+Q1_p3(i);
    P9=Qz_p3(i)+Q6_p3(i)+Q5_p3(i)+Q4_p3(i)+Q3_p3(i)+Q2_p3(i)+Q1_p3(i);
    HDFnew7212Q_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                     +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                     +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                     +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                     +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                     +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction5...
                     +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                     +((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction5);
    R1=0;
    R2=Dz_p3(i);
    R3=Dz_p3(i)+D6_p3(i);
    R4=Dz_p3(i)+D6_p3(i)+D5_p3(i);
    R5=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i);
    R6=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i);
    R7=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i)+D2_p3(i)*fraction5;
    R8=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i)+D2_p3(i)*fraction5+D1_p3(i);
    R9=Dz_p3(i)+D6_p3(i)+D5_p3(i)+D4_p3(i)+D3_p3(i)+D2_p3(i)+D1_p3(i);
    HDFnew7212D_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                     +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                     +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                     +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                     +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                     +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction5...
                     +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                     +((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction5);
end
HDF = HDFQ_p1 + HDFD_p1 + HDFQ_p2 + HDFD_p2 + HDFQ_p3 + HDFD_p3;
HDF0121 = HDF0121Q_p1 + HDF0121D_p1 + HDF0121Q_p2 + HDF0121D_p2 + HDF0121Q_p3 + HDF0121D_p3;
HDF7212 = HDF7212Q_p1 + HDF7212D_p1 + HDF7212Q_p2 + HDF7212D_p2 + HDF7212Q_p3 + HDF7212D_p3;
HDFnew0121 = HDFnew0121Q_p1 + HDFnew0121D_p1 + HDFnew0121Q_p2 + HDFnew0121D_p2 + HDFnew0121Q_p3 + HDFnew0121D_p3;
HDFnew7212 = HDFnew7212Q_p1 + HDFnew7212D_p1 + HDFnew7212Q_p2 + HDFnew7212D_p2 + HDFnew7212Q_p3 + HDFnew7212D_p3;

threeHDF=[HDF;HDF0121;HDF7212];
fiveHDF=[HDF;HDF0121;HDF7212;HDFnew0121;HDFnew7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
[minHDFvalue5,minHDFcase5]=min(fiveHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
point5=sum(minHDFcase5==4,2)+sum(minHDFcase5==5,2);
allarea3=sum(abs(HDF-minHDFvalue3));
allarea5=sum(abs(minHDFvalue3-minHDFvalue5));

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
T5bit=round(T5*bit/2);
T6bit=round(T6*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit-T5bit-T6bit;
for i=1:fsw*time%璶祘Α絏
   switch minHDFcase5(i)
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
       case 4%new01234565
           out0=repmat([0;0;0;0;0;0;0],1,Tzbit(i));
           out1=repmat(v1(:,i),1,T1bit(i));
           out2=repmat(v2(:,i),1,T2bit(i));
           out3=repmat(v3(:,i),1,T3bit(i));
           out4=repmat(v4(:,i),1,T4bit(i));
           out5=repmat(v5(:,i),1,floor(T5bit(i)*fraction5));
           out6=repmat(v6(:,i),1,T6bit(i));
           out7=repmat(v5(:,i),1,T5bit(i)-floor(T5bit(i)*fraction5));
           out01234565=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out01234565;
       case 5%new76543212
           out0=repmat([1;1;1;1;1;1;1],1,Tzbit(i));
           out1=repmat(v6(:,i),1,T6bit(i));
           out2=repmat(v5(:,i),1,T5bit(i));
           out3=repmat(v4(:,i),1,T4bit(i));
           out4=repmat(v3(:,i),1,T3bit(i));
           out5=repmat(v2(:,i),1,floor(T2bit(i)*fraction5));
           out6=repmat(v1(:,i),1,T1bit(i));
           out7=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)*fraction5));
           out76543212=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out76543212;
   end
end
out_2=circshift(out,[0 1]);
out_3=[out;out_2];
diff1=sum(abs(out(1,:)-out_2(1,:)));
diff2=sum(abs(out(2,:)-out_2(2,:)));
diff3=sum(abs(out(3,:)-out_2(3,:)));
diff4=sum(abs(out(4,:)-out_2(4,:)));
diff5=sum(abs(out(5,:)-out_2(5,:)));
diff6=sum(abs(out(6,:)-out_2(6,:)));
diff7=sum(abs(out(7,:)-out_2(7,:)));
fprintf('out1ち传Ω计=%f\n',diff1);
fprintf('out2ち传Ω计=%f\n',diff2);
fprintf('out3ち传Ω计=%f\n',diff3);
fprintf('out4ち传Ω计=%f\n',diff4);
fprintf('out5ち传Ω计=%f\n',diff5);
fprintf('out6ち传Ω计=%f\n',diff6);
fprintf('out7ち传Ω计=%f\n',diff7);

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
disperence=[point5;allarea5;THD_Van;ITHD_Van;ITHD_Vbn;ITHDavg];
% figure;%allHDF
% plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
% plot(angle*180/pi,HDF0121,'b'),hold on;
% plot(angle*180/pi,HDF7212,'g'),hold on;
% plot(angle*180/pi,minHDFvalue3,'m'),hold on;
% plot(angle*180/pi,HDFnew0121,'c'),hold on;
% plot(angle*180/pi,HDFnew7212,'y'),hold on;
% plot(angle*180/pi,minHDFvalue5,'k'),hold on;