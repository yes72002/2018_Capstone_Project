clear%11/15き穝3跋
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.52;%0.2472,0.4,0.6472
fraction3=0.9;
f = 50;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t);
s2 = amp*cos(2*pi*f*t-2*pi/5); 
s3 = amp*cos(2*pi*f*t-4*pi/5);
s4 = amp*cos(2*pi*f*t-6*pi/5); 
s5 = amp*cos(2*pi*f*t-8*pi/5); 
a = 2*pi/5;
mapping=2/5*[...
    cos(0) cos(a) cos(2*a) cos(3*a) cos(4*a);...
    sin(0) sin(a) sin(2*a) sin(3*a) sin(4*a);...
    cos(0) cos(2*a) cos(4*a) cos(6*a) cos(8*a);...
    sin(0) sin(2*a) sin(4*a) sin(6*a) sin(8*a);...
    1 1 1 1 1];
sindqxy=mapping*[s1;s2;s3;s4;s5];
sindq1x=sindqxy(1,:);
sindq1y=sindqxy(2,:);
sindq2x=sindqxy(3,:);
sindq2y=sindqxy(4,:);
angle=atan2(sindq1y,sindq1x);
for i=1:1:fsw*time
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
% plot(sindq1y,sindq1x);
% title('sindq1x,sindq1y');
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
%{
figure;
quiver(zeros(1,30),zeros(1,30),DQ(1,:),DQ(2,:))
title('dq-plane 1'); 
for i=1:30
    k=num2str(i);
    text(DQ(1,i),DQ(2,i),k);
end
figure;
quiver(zeros(1,30),zeros(1,30),DQ(3,:),DQ(4,:))
title('dq-plane 2'); 
for i=1:30
    k=num2str(i);
    text(DQ(3,i),DQ(4,i),k);
end
%}
duty=zeros(4,fsw*time);
sector=zeros(1,fsw*time);
v1=zeros(5,fsw*time);
v2=zeros(5,fsw*time);
v3=zeros(5,fsw*time);
v4=zeros(5,fsw*time);
for i=1:fsw*time%耞sector
if angle(i)>0 && angle(i)<=pi/5;%sector1
    sector(i)=1;
    duty(:,i)=[M1(:,16) M1(:,24) M1(:,25) M1(:,29);M2(:,16) M2(:,24) M2(:,25) M2(:,29)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[1;0;0;0;0];
    v2(:,i)=[1;1;0;0;0];
    v3(:,i)=[1;1;0;0;1];
    v4(:,i)=[1;1;1;0;1];
end
if angle(i)>pi/5 && angle(i)<=2*pi/5;%sector2
    sector(i)=2;
    duty(:,i)=[M1(:,8) M1(:,24) M1(:,28) M1(:,29);M2(:,8) M2(:,24) M2(:,28) M2(:,29)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;1;0;0;0];
    v2(:,i)=[1;1;0;0;0];
    v3(:,i)=[1;1;1;0;0];
    v4(:,i)=[1;1;1;0;1];
end   
if angle(i)>2*pi/5 && angle(i)<=3*pi/5;%sector3
    sector(i)=3;
    duty(:,i)=[M1(:,8) M1(:,12) M1(:,28) M1(:,30);M2(:,8) M2(:,12) M2(:,28) M2(:,30)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;1;0;0;0];
    v2(:,i)=[0;1;1;0;0];
    v3(:,i)=[1;1;1;0;0];
    v4(:,i)=[1;1;1;1;0];
end
if angle(i)>3*pi/5 && angle(i)<=4*pi/5;%sector4
    sector(i)=4;
    duty(:,i)=[M1(:,4) M1(:,12) M1(:,14) M1(:,30);M2(:,4) M2(:,12) M2(:,14) M2(:,30)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;0;1;0;0];
    v2(:,i)=[0;1;1;0;0];
    v3(:,i)=[0;1;1;1;0];
    v4(:,i)=[1;1;1;1;0];
end
if angle(i)>4*pi/5 && angle(i)<=pi;%sector5
    sector(i)=5;
    duty(:,i)=[M1(:,4) M1(:,6) M1(:,14) M1(:,15);M2(:,4) M2(:,6) M2(:,14) M2(:,15)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;0;1;0;0];
    v2(:,i)=[0;0;1;1;0];
    v3(:,i)=[0;1;1;1;0];
    v4(:,i)=[0;1;1;1;1];
end 
if angle(i)>pi && angle(i)<=6*pi/5;%sector6
    sector(i)=6;
    duty(:,i)=[M1(:,2) M1(:,6) M1(:,7) M1(:,15);M2(:,2) M2(:,6) M2(:,7) M2(:,15)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;0;0;1;0];
    v2(:,i)=[0;0;1;1;0];
    v3(:,i)=[0;0;1;1;1];
    v4(:,i)=[0;1;1;1;1];
end
if angle(i)>6*pi/5 && angle(i)<=7*pi/5;%sector7
    sector(i)=7;
    duty(:,i)=[M1(:,2) M1(:,3) M1(:,7) M1(:,23);M2(:,2) M2(:,3) M2(:,7) M2(:,23)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;0;0;1;0];
    v2(:,i)=[0;0;0;1;1];
    v3(:,i)=[0;0;1;1;1];
    v4(:,i)=[1;0;1;1;1];
end
if angle(i)>7*pi/5 && angle(i)<=8*pi/5;%sector8
    sector(i)=8;
    duty(:,i)=[M1(:,1) M1(:,3) M1(:,19) M1(:,23);M2(:,1) M2(:,3) M2(:,19) M2(:,23)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;0;0;0;1];
    v2(:,i)=[0;0;0;1;1];
    v3(:,i)=[1;0;0;1;1];
    v4(:,i)=[1;0;1;1;1];
end
if angle(i)>8*pi/5 && angle(i)<=9*pi/5;%sector9
    sector(i)=9;
    duty(:,i)=[M1(:,1) M1(:,17) M1(:,19) M1(:,27);M2(:,1) M2(:,17) M2(:,19) M2(:,27)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[0;0;0;0;1];
    v2(:,i)=[1;0;0;0;1];
    v3(:,i)=[1;0;0;1;1];
    v4(:,i)=[1;1;0;1;1];
end
if angle(i)>9*pi/5 && angle(i)<=2*pi;%sector10
    sector(i)=10;
    duty(:,i)=[M1(:,16) M1(:,17) M1(:,25) M1(:,27);M2(:,16) M2(:,17) M2(:,25) M2(:,27)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    v1(:,i)=[1;0;0;0;0];
    v2(:,i)=[1;0;0;0;1];
    v3(:,i)=[1;1;0;0;1];
    v4(:,i)=[1;1;0;1;1];
end
end
a1=duty(1,:);
a2=duty(2,:);
a3=duty(3,:);
a4=duty(4,:);
az=1-a1-a2-a3-a4;
a0=az/2;
a7=a0;
% k=[sector;a0;a1;a2;a3;a4;a7];

T1=a1;%礚蛤Τ弘非
T2=a2;
T3=a3;
T4=a4;
Tz=az;
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(i),2)==1%1,sector碞琌计0,sector碞琌案计
        Q1(i)=((M1(1,16))*cos(angle(i)-(pi/5)*(sector(i)-1))-amp)*T1(i);
        Q2(i)=((M1(1,25))*cos((pi/5)*sector(i)-angle(i))-amp)*T2(i);
        Q3(i)=((M1(1,25))*cos(angle(i)-(pi/5)*(sector(i)-1))-amp)*T3(i);
        Q4(i)=((M1(1,16))*cos((pi/5)*sector(i)-angle(i))-amp)*T4(i);
        Qz(i)=-amp*Tz(i);
        D1(i)=((M1(1,16))*sin(angle(i)-(pi/5)*(sector(i)-1)))*T1(i);
        D2(i)=((-M1(1,25))*sin((pi/5)*sector(i)-angle(i)))*T2(i);
        D3(i)=((M1(1,25))*sin(angle(i)-(pi/5)*(sector(i)-1)))*T3(i);
        D4(i)=((-M1(1,16))*sin((pi/5)*sector(i)-angle(i)))*T4(i);
        Dz(i)=0;
    else
        Q1(i)=((M1(1,16))*cos((pi/5)*sector(i)-angle(i))-amp)*T1(i);
        Q2(i)=((M1(1,25))*cos(angle(i)-(pi/5)*(sector(i)-1))-amp)*T2(i);
        Q3(i)=((M1(1,25))*cos((pi/5)*sector(i)-angle(i))-amp)*T3(i);
        Q4(i)=((M1(1,16))*cos(angle(i)-(pi/5)*(sector(i)-1))-amp)*T4(i);
        Qz(i)=-amp*Tz(i); 
        D1(i)=((M1(1,16))*sin((pi/5)*sector(i)-angle(i)))*T1(i);
        D2(i)=((-M1(1,25))*sin(angle(i)-(pi/5)*(sector(i)-1)))*T2(i);
        D3(i)=((M1(1,25))*sin((pi/5)*sector(i)-angle(i)))*T3(i);
        D4(i)=((-M1(1,16))*sin(angle(i)-(pi/5)*(sector(i)-1)))*T4(i);
        Dz(i)=0;
    end
end
QQQ=Qz+Q1+Q2+Q3+Q4;
DDD=Dz+D1+D2+D3+D4;
for i=1:fsw*time
    P1 = 0;%HDF012347
    P2 = Qz(i)/2;
    P3 = Qz(i)/2+Q1(i);
    P4 = Qz(i)/2+Q1(i)+Q2(i);
    P5 = Qz(i)/2+Q1(i)+Q2(i)+Q3(i);
    P6 = Qz(i)/2+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    P7 = Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    HDFQ(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*Tz(i)/2;
    R1 = 0;
    R2 = Dz(i)/2;
    R3 = Dz(i)/2+D1(i);
    R4 = Dz(i)/2+D1(i)+D2(i);
    R5 = Dz(i)/2+D1(i)+D2(i)+D3(i);
    R6 = Dz(i)/2+D1(i)+D2(i)+D3(i)+D4(i);
    R7 = Dz(i)+D1(i)+D2(i)+D3(i)+D4(i);
    HDFD(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*Tz(i)/2;
    P1 = 0;%HDF012343
    P2 = Qz(i);
    P3 = Qz(i)+Q1(i);
    P4 = Qz(i)+Q1(i)+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i)+Q3(i)*fraction3;
    P6 = Qz(i)+Q1(i)+Q2(i)+Q3(i)*fraction3+Q4(i);
    P7 = Qz(i)+Q1(i)+Q2(i)+Q3(i)+Q4(i);
    HDF0121Q(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(i)*(1-fraction3);
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D1(i);
    R4 = Dz(i)+D1(i)+D2(i);
    R5 = Dz(i)+D1(i)+D2(i)+D3(i)*fraction3;
    R6 = Dz(i)+D1(i)+D2(i)+D3(i)*fraction3+D4(i);
    R7 = Dz(i)+D1(i)+D2(i)+D3(i)+D4(i);
    HDF0121D(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Qz(i);
    P3 = Qz(i)+Q4(i);
    P4 = Qz(i)+Q4(i)+Q3(i);
    P5 = Qz(i)+Q4(i)+Q3(i)+Q2(i)*fraction3;
    P6 = Qz(i)+Q4(i)+Q3(i)+Q2(i)*fraction3+Q1(i);
    P7 = Qz(i)+Q4(i)+Q3(i)+Q2(i)+Q1(i);
    HDF7212Q(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D4(i);
    R4 = Dz(i)+D4(i)+D3(i);
    R5 = Dz(i)+D4(i)+D3(i)+D2(i)*fraction3;
    R6 = Dz(i)+D4(i)+D3(i)+D2(i)*fraction3+D1(i);
    R7 = Dz(i)+D4(i)+D3(i)+D2(i)+D1(i);
    HDF7212D(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*(1-fraction3);
end
HDF = HDFQ + HDFD;
HDF0121 = HDF0121Q + HDF0121D;
HDF7212 = HDF7212Q + HDF7212D;

threeHDF=[HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
allarea3=sum(abs(HDF-minHDFvalue3));

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit;
for i=1:fsw*time%璶祘Α絏
   switch minHDFcase3(i)
       case 1%012347
           out0=repmat([0;0;0;0;0],1,floor(Tzbit(i)/2));
           out1=repmat(v1(:,i),1,T1bit(i));
           out2=repmat(v2(:,i),1,T2bit(i));
           out3=repmat(v3(:,i),1,T3bit(i));
           out4=repmat(v4(:,i),1,T4bit(i));
           out7=repmat([1;1;1;1;1],1,Tzbit(i)-floor(Tzbit(i)/2));
           out012347=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out012347;
       case 2%012343
           out0=repmat([0;0;0;0;0],1,Tzbit(i));
           out1=repmat(v1(:,i),1,T1bit(i));
           out2=repmat(v2(:,i),1,T2bit(i));
           out3=repmat(v3(:,i),1,floor(T3bit(i)*fraction3));
           out4=repmat(v4(:,i),1,T4bit(i));
           out7=repmat(v3(:,i),1,T3bit(i)-floor(T3bit(i)*fraction3));
           out012343=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out012343;
       case 3%743212
           out0=repmat([1;1;1;1;1],1,Tzbit(i));
           out1=repmat(v4(:,i),1,T4bit(i));
           out2=repmat(v3(:,i),1,T3bit(i));
           out3=repmat(v2(:,i),1,floor(T2bit(i)*fraction3));
           out4=repmat(v1(:,i),1,T1bit(i));
           out7=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)*fraction3));
           out743212=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=out743212;
   end
end
out_2=circshift(out,[0 1]);
out3=[out;out_2];
diff1=sum(abs(out(1,:)-out_2(1,:)));
diff2=sum(abs(out(2,:)-out_2(2,:)));
diff3=sum(abs(out(3,:)-out_2(3,:)));
diff4=sum(abs(out(4,:)-out_2(4,:)));
diff5=sum(abs(out(5,:)-out_2(5,:)));
fprintf('out1ち传Ω计=%f\n',diff1);
fprintf('out2ち传Ω计=%f\n',diff2);
fprintf('out3ち传Ω计=%f\n',diff3);
fprintf('out4ち传Ω计=%f\n',diff4);
fprintf('out5ち传Ω计=%f\n',diff5);

N = time*fsw*bit;
Van=[4/5 -1/5 -1/5 -1/5 -1/5]*out;
Vbn=[-1/5 4/5 -1/5 -1/5 -1/5]*out;
Vcn=[-1/5 -1/5 4/5 -1/5 -1/5]*out;
Vdn=[-1/5 -1/5 -1/5 4/5 -1/5]*out;
Ven=[-1/5 -1/5 -1/5 -1/5 4/5]*out;
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
k=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:い丁
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(繵^2/膀繵^2)
THD_Van=100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
THD_Vdn=100*((sum(mag_Vdn_1.^2)-(mag_Vdn_1(1,1).^2))/(mag_Vdn_1(1,1).^2)).^(1/2);
THD_Ven=100*((sum(mag_Ven_1.^2)-(mag_Ven_1(1,1).^2))/(mag_Ven_1(1,1).^2)).^(1/2);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n',THD_Vcn);
fprintf('VTHD4=%f\n',THD_Vdn);
fprintf('VTHD5=%f\n',THD_Ven);
%ITHD=繵瞯埃セō
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%ソ搭埃そ畉+1=[(N/2-4)-6]/5+1
for i=1:length(mag_Van_1)
    mag_ithd_out1(i)=mag_Van_1(i)/(10*i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(10*i); 
    mag_ithd_out3(i)=mag_Vcn_1(i)/(10*i);  
    mag_ithd_out4(i)=mag_Vdn_1(i)/(10*i);  
    mag_ithd_out5(i)=mag_Ven_1(i)/(10*i);  
end
ITHD_Van=100*((sum(mag_ithd_out1.^2)-(mag_ithd_out1(1,1).^2))/(mag_ithd_out1(1,1).^2)).^(1/2);
ITHD_Vbn=100*((sum(mag_ithd_out2.^2)-(mag_ithd_out2(1,1).^2))/(mag_ithd_out2(1,1).^2)).^(1/2);
ITHD_Vcn=100*((sum(mag_ithd_out3.^2)-(mag_ithd_out3(1,1).^2))/(mag_ithd_out3(1,1).^2)).^(1/2);
ITHD_Vdn=100*((sum(mag_ithd_out4.^2)-(mag_ithd_out4(1,1).^2))/(mag_ithd_out4(1,1).^2)).^(1/2);
ITHD_Ven=100*((sum(mag_ithd_out5.^2)-(mag_ithd_out5(1,1).^2))/(mag_ithd_out5(1,1).^2)).^(1/2);
fprintf('ITHD1=%f\n',ITHD_Van);
fprintf('ITHD2=%f\n',ITHD_Vbn);
fprintf('ITHD3=%f\n',ITHD_Vcn);
fprintf('ITHD4=%f\n',ITHD_Vdn);
fprintf('ITHD5=%f\n',ITHD_Ven);
%赣衡,︽常衡ㄓ,钡ㄓ琌礶瓜
n2=2:1:N/2;%繷(繷材Τ粇箂)
mag_Van_2(1,n2-1)=mag_Van(1,n2);%5:0.499,10:0.00007,15...
mag_Vbn_2(1,n2-1)=mag_Vbn(1,n2); 
mag_Vcn_2(1,n2-1)=mag_Vcn(1,n2); 
mag_Vdn_2(1,n2-1)=mag_Vdn(1,n2); 
mag_Ven_2(1,n2-1)=mag_Ven(1,n2); 
mag_Van_3=zeros(1,(N/2-1)*10);
mag_Vbn_3=zeros(1,(N/2-1)*10);
mag_Vcn_3=zeros(1,(N/2-1)*10);
mag_Vdn_3=zeros(1,(N/2-1)*10);
mag_Ven_3=zeros(1,(N/2-1)*10);
for k=1:1:N/2-1;%5->50
mag_Van_3(1,(bit*fsw/N)*k)=mag_Van_2(1,k);
mag_Vbn_3(1,(bit*fsw/N)*k)=mag_Vbn_2(1,k);
mag_Vcn_3(1,(bit*fsw/N)*k)=mag_Vcn_2(1,k);
mag_Vdn_3(1,(bit*fsw/N)*k)=mag_Vdn_2(1,k);
mag_Ven_3(1,(bit*fsw/N)*k)=mag_Ven_2(1,k);
end
% figure===================================================================
% semilogx(mag_Van_3,'r'),title('mag Van 3');
end_time=clock;
execution_time=end_time-start_time;
ITHDavg=(ITHD_Van+ITHD_Vbn+ITHD_Vcn+ITHD_Vdn+ITHD_Ven)/5;
disperence=[point3;allarea3;THD_Van;ITHD_Van;ITHDavg];
figure;%allHDF
degree=angle*180/pi;
plot(degree,HDF,'r'),title('allHDF'),hold on;
plot(degree,HDF0121,'b'),hold on;
plot(degree,HDF7212,'g'),hold on;
plot(degree,minHDFvalue3,'m'),hold on;