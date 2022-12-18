clear%10/18新3區
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.5;
fraction3=0.5;
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
for i=1:1:fsw*time%HDF0127
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
    P1 = 0;%HDF0121
    P2 = Qz(i);
    P3 = Qz(i)+Q1(i)*fraction3;
    P4 = Qz(i)+Q1(i)*fraction3+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i);
    HDF0121Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)*fraction3;
    HDF0121Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDF0121Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T1(i)*(1-fraction3);    
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D1(i)*fraction3;
    R4 = Dz(i)+D1(i)*fraction3+D2(i);
    R5 = Dz(i)+D1(i)+D2(i);
    HDF0121D1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF0121D2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)*fraction3;
    HDF0121D3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDF0121D4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T1(i)*(1-fraction3);
    P1 = 0;%HDF7212
    P2 = Qz(i);
    P3 = Qz(i)+Q2(i)*fraction3;
    P4 = Qz(i)+Q2(i)*fraction3+Q1(i);
    P5 = Qz(i)+Q2(i)+Q1(i);
    HDF7212Q1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF7212Q2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T2(i)*fraction3;
    HDF7212Q3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T1(i);
    HDF7212Q4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D2(i)*fraction3;
    R4 = Dz(i)+D2(i)*fraction3+D1(i);
    R5 = Dz(i)+D2(i)+D1(i);
    HDF7212D1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF7212D2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T2(i)*fraction3;
    HDF7212D3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T1(i);
    HDF7212D4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*(1-fraction3);
end
HDFQ=HDFQ1+HDFQ2+HDFQ3+HDFQ4;
HDFD=HDFD1+HDFD2+HDFD3+HDFD4;
HDF=HDFQ+HDFD;
HDF0121Q=HDF0121Q1+HDF0121Q2+HDF0121Q3+HDF0121Q4;
HDF0121D=HDF0121D1+HDF0121D2+HDF0121D3+HDF0121D4;
HDF0121=HDF0121Q+HDF0121D;
HDF7212Q=HDF7212Q1+HDF7212Q2+HDF7212Q3+HDF7212Q4;
HDF7212D=HDF7212D1+HDF7212D2+HDF7212D3+HDF7212D4;
HDF7212=HDF7212Q+HDF7212D;

threeHDF=[HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
allarea3=sum(abs(HDF-minHDFvalue3));

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
           out001=repmat(v1(:,i),1,floor(T1bit(i)*fraction3));
           out011=repmat(v2(:,i),1,T2bit(i));
           out111=repmat(v1(:,i),1,T1bit(i)-floor(T1bit(i)*fraction3));
           out0121=[out000 out001 out011 out111 out111 out011 out001 out000];
           out(:,bit*(i-1)+1:bit*i)=out0121;
       case 3%7212
           out000=repmat([1;1;1],1,Tzbit(i));
           out001=repmat(v2(:,i),1,floor(T2bit(i)*fraction3));
           out011=repmat(v1(:,i),1,T1bit(i));
           out111=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)*fraction3));
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

N = time*fsw*bit;
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
% figure;%=================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
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
    mag_ithd_out1(i)=mag_Van_1(i)/(i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(i); 
    mag_ithd_out3(i)=mag_Vcn_1(i)/(i);  
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
% figure;
% semilogx(mag_Van_3,'r'),title('mag Van 3');
end_time=clock;
execution_time=end_time-start_time;
ITHDavg=(ITHD_Van+ITHD_Vbn+ITHD_Vcn)/3;
disperence=[point3;allarea3;THD_Van;ITHD_Van;ITHDavg];
% figure;%allHDF
% plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
% plot(angle*180/pi,HDF0121,'b'),hold on;
% plot(angle*180/pi,HDF7212,'g'),hold on;
% plot(angle*180/pi,minHDFvalue3,'m'),hold on;