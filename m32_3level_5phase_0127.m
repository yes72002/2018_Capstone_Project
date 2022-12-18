clear%3/3三相
close all
clc
start_time = clock;
time  =  0.1;
fsw  =  5000;
bit  =  400;
amp  =  0.5;%0.2472,0.4,0.6472
f  =  50;
n  =  fsw*time;
t  =  0:1/fsw:time-1/fsw; 
s1  =  amp*cos(2*pi*f*t);
s2  =  amp*cos(2*pi*f*t-2*pi/5); 
s3  =  amp*cos(2*pi*f*t-4*pi/5);
s4  =  amp*cos(2*pi*f*t-6*pi/5); 
s5  =  amp*cos(2*pi*f*t-8*pi/5);
allsin = [s1;s2;s3;s4;s5];
duty = allsin-floor(allsin);
dutybit = duty*bit;
updown = floor(allsin);
for i = 1:fsw*time%主要的大程式碼
    ud1 = updown(1,i);
    ud2 = updown(2,i);
    ud3 = updown(3,i);
    ud4 = updown(4,i);
    ud5 = updown(5,i);
    o1 = round(dutybit(1,i));
    o2 = round(dutybit(2,i));
    o3 = round(dutybit(3,i));
    o4 = round(dutybit(4,i));
    o5 = round(dutybit(5,i));
    z1 = (bit-o1)/2;
    z2 = (bit-o2)/2;
    z3 = (bit-o3)/2;
    z4 = (bit-o4)/2;
    z5 = (bit-o5)/2;
    out1(:,bit*(i-1)+1:bit*i)=[zeros(1,floor(z1))+ud1 ones(1,o1)+ud1 zeros(1,ceil(z1))+ud1];
    out2(:,bit*(i-1)+1:bit*i)=[zeros(1,floor(z2))+ud2 ones(1,o2)+ud2 zeros(1,ceil(z2))+ud2];
    out3(:,bit*(i-1)+1:bit*i)=[zeros(1,floor(z3))+ud3 ones(1,o3)+ud3 zeros(1,ceil(z3))+ud3];
    out4(:,bit*(i-1)+1:bit*i)=[zeros(1,floor(z4))+ud4 ones(1,o4)+ud4 zeros(1,ceil(z4))+ud4];
    out5(:,bit*(i-1)+1:bit*i)=[zeros(1,floor(z5))+ud5 ones(1,o5)+ud5 zeros(1,ceil(z5))+ud5];
end
out=[out1;out2;out3;out4;out5];
figure;
subplot(6,1,1),plot(t,s1,t,s2,t,s3,t,s4,t,s5);
subplot(6,1,2),plot(out(1,:));
subplot(6,1,3),plot(out(2,:));
subplot(6,1,4),plot(out(3,:));
subplot(6,1,5),plot(out(4,:));
subplot(6,1,6),plot(out(5,:));

N = time*fsw*bit;
Van=[4/5 -1/5 -1/5 -1/5 -1/5]*out;
Vbn=[-1/5 4/5 -1/5 -1/5 -1/5]*out;
Vcn=[-1/5 -1/5 4/5 -1/5 -1/5]*out;
Vdn=[-1/5 -1/5 -1/5 4/5 -1/5]*out;
Ven=[-1/5 -1/5 -1/5 -1/5 4/5]*out;
fft_Van=fft(Van,N)/N;%6:0.499且兩端都有值
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
%i=6:5:到中間
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(倍頻^2/基頻^2)
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
%ITHD=頻率除以本身值
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%末減首除以公差+1=[(N/2-4)-6]/5+1
for i=1:length(mag_Van_1)
    mag_ithd_out1(i)=mag_Van_1(i)/(i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(i); 
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
%該算的,此行以上都算出來了,接下來是畫圖------------------------------------
n2=2:1:N/2;%去頭(頭第一個值有誤為零)
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
