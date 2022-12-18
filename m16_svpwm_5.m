clear%8/1,8/2,8/4き1跋
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.52;%0.2472,0.4,0.6472
f = 50;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t);
s2 = amp*cos(2*pi*f*t-2*pi/5); 
s3 = amp*cos(2*pi*f*t-4*pi/5);
s4 = amp*cos(2*pi*f*t-6*pi/5); 
s5 = amp*cos(2*pi*f*t-8*pi/5); 
a=2*pi/5;
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
figure;
plot(sindq1y,sindq1x);
hold on;
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

duty=zeros(4,fsw*time);
sector=zeros(1,fsw*time);
a0=zeros(1,fsw*time);
a7=zeros(1,fsw*time);
a1=zeros(1,fsw*time);
a2=zeros(1,fsw*time);
a3=zeros(1,fsw*time);
a4=zeros(1,fsw*time);
v1=zeros(5,fsw*time);
v2=zeros(5,fsw*time);
v3=zeros(5,fsw*time);
v4=zeros(5,fsw*time);
for i=1:fsw*time%耞sector%jim
if angle(i)>0 && angle(i)<=pi/5;%sector1
    sector(i)=1;
    duty=[M1(:,16) M1(:,24) M1(:,25) M1(:,29);M2(:,16) M2(:,24) M2(:,25) M2(:,29)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[1;0;0;0;0];
    v2(:,i)=[1;1;0;0;0];
    v3(:,i)=[1;1;0;0;1];
    v4(:,i)=[1;1;1;0;1];
end
if angle(i)>pi/5 && angle(i)<=2*pi/5;%sector2
    sector(i)=2;
    duty=[M1(:,8) M1(:,24) M1(:,28) M1(:,29);M2(:,8) M2(:,24) M2(:,28) M2(:,29)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;1;0;0;0];
    v2(:,i)=[1;1;0;0;0];
    v3(:,i)=[1;1;1;0;0];
    v4(:,i)=[1;1;1;0;1];
end   
if angle(i)>2*pi/5 && angle(i)<=3*pi/5;%sector3
    sector(i)=3;
    duty=[M1(:,8) M1(:,12) M1(:,28) M1(:,30);M2(:,8) M2(:,12) M2(:,28) M2(:,30)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;1;0;0;0];
    v2(:,i)=[0;1;1;0;0];
    v3(:,i)=[1;1;1;0;0];
    v4(:,i)=[1;1;1;1;0];
end
if angle(i)>3*pi/5 && angle(i)<=4*pi/5;%sector4
    sector(i)=4;
    duty=[M1(:,4) M1(:,12) M1(:,14) M1(:,30);M2(:,4) M2(:,12) M2(:,14) M2(:,30)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;0;1;0;0];
    v2(:,i)=[0;1;1;0;0];
    v3(:,i)=[0;1;1;1;0];
    v4(:,i)=[1;1;1;1;0];
end
if angle(i)>4*pi/5 && angle(i)<=pi;%sector5
    sector(i)=5;
    duty=[M1(:,4) M1(:,6) M1(:,14) M1(:,15);M2(:,4) M2(:,6) M2(:,14) M2(:,15)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;0;1;0;0];
    v2(:,i)=[0;0;1;1;0];
    v3(:,i)=[0;1;1;1;0];
    v4(:,i)=[0;1;1;1;1];
end 
if angle(i)>pi && angle(i)<=6*pi/5;%sector6
    sector(i)=6;
    duty=[M1(:,2) M1(:,6) M1(:,7) M1(:,15);M2(:,2) M2(:,6) M2(:,7) M2(:,15)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;0;0;1;0];
    v2(:,i)=[0;0;1;1;0];
    v3(:,i)=[0;0;1;1;1];
    v4(:,i)=[0;1;1;1;1];
end
if angle(i)>6*pi/5 && angle(i)<=7*pi/5;%sector7
    sector(i)=7;
    duty=[M1(:,2) M1(:,3) M1(:,7) M1(:,23);M2(:,2) M2(:,3) M2(:,7) M2(:,23)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;0;0;1;0];
    v2(:,i)=[0;0;0;1;1];
    v3(:,i)=[0;0;1;1;1];
    v4(:,i)=[1;0;1;1;1];
end
if angle(i)>7*pi/5 && angle(i)<=8*pi/5;%sector8
    sector(i)=8;
    duty=[M1(:,1) M1(:,3) M1(:,19) M1(:,23);M2(:,1) M2(:,3) M2(:,19) M2(:,23)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;0;0;0;1];
    v2(:,i)=[0;0;0;1;1];
    v3(:,i)=[1;0;0;1;1];
    v4(:,i)=[1;0;1;1;1];
end
if angle(i)>8*pi/5 && angle(i)<=9*pi/5;%sector9
    sector(i)=9;
    duty=[M1(:,1) M1(:,17) M1(:,19) M1(:,27);M2(:,1) M2(:,17) M2(:,19) M2(:,27)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[0;0;0;0;1];
    v2(:,i)=[1;0;0;0;1];
    v3(:,i)=[1;0;0;1;1];
    v4(:,i)=[1;1;0;1;1];
end
if angle(i)>9*pi/5 && angle(i)<=2*pi;%sector10
    sector(i)=10;
    duty=[M1(:,16) M1(:,17) M1(:,25) M1(:,27);M2(:,16) M2(:,17) M2(:,25) M2(:,27)]\[sindq1x(i);sindq1y(i);sindq2x(i);sindq2y(i)];
    a1(i)=duty(1);
    a2(i)=duty(2);
    a3(i)=duty(3);
    a4(i)=duty(4);
    a0(i)=(1-a1(i)-a2(i)-a3(i)-a4(i))/2;
    a7(i)=a0(i);
    v1(:,i)=[1;0;0;0;0];
    v2(:,i)=[1;0;0;0;1];
    v3(:,i)=[1;1;0;0;1];
    v4(:,i)=[1;1;0;1;1];
end
end
k=[sector;a0;a1;a2;a3;a4;a7];

T1=a1;%礚蛤Τ弘非
T2=a2;
T3=a3;
T4=a4;
T0=a0;
T7=a7;
az=a0+a7;
Tz=az;
T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit;
TTT=T1bit+T2bit+T3bit+T4bit+Tzbit;
TTTTT=[T1bit;T2bit;T3bit;T4bit;Tzbit;TTT];
out=zeros(5,fsw*time);
for i=1:fsw*time%璶祘Α絏
    out0=repmat([0;0;0;0;0],1,floor(Tzbit(i)/2));
    out1=repmat(v1(:,i),1,T1bit(i));
    out2=repmat(v2(:,i),1,T2bit(i));
    out3=repmat(v3(:,i),1,T3bit(i));
    out4=repmat(v4(:,i),1,T4bit(i));
    out7=repmat([1;1;1;1;1],1,ceil(Tzbit(i)/2));
    out012347=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
    out(:,bit*(i-1)+1:bit*i)=out012347;
end
out_2(1,:)=out(1,2:end);
out_2(2,:)=out(2,2:end);
out_2(3,:)=out(3,2:end);
out_2(4,:)=out(4,2:end);
out_2(5,:)=out(5,2:end);
out_2(1,fsw*time*bit)=0;
out_2(2,fsw*time*bit)=0;
out_2(3,fsw*time*bit)=0;
out_2(4,fsw*time*bit)=0;
out_2(5,fsw*time*bit)=0;
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
fprintf('out5ち传Ω计=%f\n\n',diff5);

figure;
subplot(6,1,1),plot(t,s1,t,s2,t,s3,t,s4,t,s5);
title('き块');
subplot(6,1,2),plot(out(1,:));
title('out1块');
subplot(6,1,3),plot(out(2,:));
title('out2块');
subplot(6,1,4),plot(out(3,:));
title('out3块');
subplot(6,1,5),plot(out(4,:));
title('out4块');
subplot(6,1,6),plot(out(5,:));
title('out5块');

N = time*fsw*bit;
Van=[4/5 -1/5 -1/5 -1/5 -1/5]*out;
Vbn=[-1/5 4/5 -1/5 -1/5 -1/5]*out;
Vcn=[-1/5 -1/5 4/5 -1/5 -1/5]*out;
Vdn=[-1/5 -1/5 -1/5 4/5 -1/5]*out;
Ven=[-1/5 -1/5 -1/5 -1/5 4/5]*out;
fft_Van=fft(Van,N)/N;
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
%k=50/10+1:50/10:100000
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
figure;
semilogx(mag_Van_1,'r');
title('mag_Van_1');

THD_Van=100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
THD_Vdn=100*((sum(mag_Vdn_1.^2)-(mag_Vdn_1(1,1).^2))/(mag_Vdn_1(1,1).^2)).^(1/2);
THD_Ven=100*((sum(mag_Ven_1.^2)-(mag_Ven_1(1,1).^2))/(mag_Ven_1(1,1).^2)).^(1/2);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n',THD_Vcn);
fprintf('VTHD4=%f\n',THD_Vdn);
fprintf('VTHD5=%f\n\n',THD_Ven);

for i=1:fsw*bit*0.02/2-1
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
fprintf('ITHD5=%f\n\n',ITHD_Ven);

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
title('mag_Van_3');
% hold on;
% plot(mag_Vcn_3,'k');
end_time=clock;
execution_time=end_time-start_time;
ITHDavg=(ITHD_Van+ITHD_Vbn+ITHD_Vcn+ITHD_Vdn+ITHD_Ven)/5;