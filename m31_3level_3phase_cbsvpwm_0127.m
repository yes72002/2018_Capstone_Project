clear%3/29三相
close all
clc
start_time = clock;
time  =  0.02;
fsw  =  3e4;
bit  =  400;
amp  =  1.14;%0.2472,0.4,0.6472
f  =  50;
n  =  fsw*time;
t  =  0:1/fsw:time-1/fsw; 
s1  =  amp*cos(2*pi*f*t);
s2  =  amp*cos(2*pi*f*t-2*pi/3); 
s3  =  amp*cos(2*pi*f*t-4*pi/3);
allsin = [s1;s2;s3];
sinmax = max(allsin);
sinmin = min(allsin);
w = (-sinmin-sinmax)/2;
allsin = allsin+[w;w;w];
duty = allsin-floor(allsin);
[dutysort,dutycase] = sort(duty);
duty3rd(1,:) = dutysort(1,:);
duty2nd(1,:) = dutysort(2,:);
duty1st(1,:) = dutysort(3,:);
Duty1 = duty1st-duty2nd;
Duty2 = duty2nd-duty3rd;
Dutyz = 1-Duty1-Duty2;
for i = 1:fsw*time
    v(dutycase(1,i),:) = [0 0 0 1];
    v(dutycase(2,i),:) = [0 0 1 1];
    v(dutycase(3,i),:) = [0 1 1 1];
    v1(:,i) = v(:,2);
    v2(:,i) = v(:,3);
end
T1 = Duty1;%無限跟有限精準度
T2 = Duty2;
Tz = Dutyz;

T1bit = round(T1*bit/2);
T2bit = round(T2*bit/2);
Tzbit = bit/2-T1bit-T2bit;
updown = floor(allsin);
for i = 1:fsw*time%主要的大程式碼
    ud=updown(:,i);
    out000 = repmat([0;0;0]+ud,1,floor(Tzbit(i)/2));
    out001 = repmat(v1(:,i)+ud,1,T1bit(i));
    out011 = repmat(v2(:,i)+ud,1,T2bit(i));
    out111 = repmat([1;1;1]+ud,1,Tzbit(i)-floor(Tzbit(i)/2));
    out0127 = [out000 out001 out011 out111 out111 out011 out001 out000];
    out(:,bit*(i-1)+1:bit*i) = out0127;
end
out_2 = circshift(out,[0 1]);
out3 = [out;out_2];
diff1 = sum(abs(out(1,:)-out_2(1,:)));
diff2 = sum(abs(out(2,:)-out_2(2,:)));
diff3 = sum(abs(out(3,:)-out_2(3,:)));
fprintf('out1切換次數 = %f\n',diff1);
fprintf('out2切換次數 = %f\n',diff2);
fprintf('out3切換次數 = %f\n',diff3);

N  =  time*fsw*bit;
Van = [ 2/3 -1/3 -1/3]*out;
Vbn = [-1/3  2/3 -1/3]*out;
Vcn = [-1/3 -1/3  2/3]*out;
fft_Van = fft(Van,N)/N;%6:0.499且兩端都有值
mag_Van = abs(fft_Van)*2;
fft_Vbn = fft(Vbn,N)/N;
mag_Vbn = abs(fft_Vbn)*2;
fft_Vcn = fft(Vcn,N)/N;
mag_Vcn = abs(fft_Vcn)*2;
i = f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i = (50/10)+1:(50/10):100000
%i = 6:5:到中間
mag_Van_1 = mag_Van(1,i);
mag_Vbn_1 = mag_Vbn(1,i);
mag_Vcn_1 = mag_Vcn(1,i);
% figure;%===========================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD = 100*sqrt(倍頻^2/基頻^2)
THD_Van = 100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn = 100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn = 100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
fprintf('VTHD1 = %f\n',THD_Van);
fprintf('VTHD2 = %f\n',THD_Vbn);
fprintf('VTHD3 = %f\n',THD_Vcn);
%ITHD = 頻率除以本身值
mag_ithd_out1 = zeros(1,N/5/2-1);
mag_ithd_out2 = zeros(1,N/5/2-1);
mag_ithd_out3 = zeros(1,N/5/2-1);
for i = 1:fsw*bit*0.02/2-1
    mag_ithd_out1(i) = mag_Van_1(i)/(i); 
    mag_ithd_out2(i) = mag_Vbn_1(i)/(i); 
    mag_ithd_out3(i) = mag_Vcn_1(i)/(i);  
end
ITHD_Van = 100*((sum(mag_ithd_out1.^2)-(mag_ithd_out1(1,1).^2))/(mag_ithd_out1(1,1).^2)).^(1/2);
ITHD_Vbn = 100*((sum(mag_ithd_out2.^2)-(mag_ithd_out2(1,1).^2))/(mag_ithd_out2(1,1).^2)).^(1/2);
ITHD_Vcn = 100*((sum(mag_ithd_out3.^2)-(mag_ithd_out3(1,1).^2))/(mag_ithd_out3(1,1).^2)).^(1/2);
fprintf('ITHD1 = %f\n',ITHD_Van);
fprintf('ITHD2 = %f\n',ITHD_Vbn);
fprintf('ITHD3 = %f\n',ITHD_Vcn);
%該算的,此行以上都算出來了,接下來是畫圖
n2 = 2:1:N/2;%去頭(頭第一個值有誤為零)
mag_Van_2(1,n2-1) = mag_Van(1,n2);%5:0.499,10:0.00007,15...
mag_Vbn_2(1,n2-1) = mag_Vbn(1,n2);
mag_Vcn_2(1,n2-1) = mag_Vcn(1,n2); 
mag_Van_3 = zeros(1,(N/2-1)*10);
mag_Vbn_3 = zeros(1,(N/2-1)*10);
mag_Vcn_3 = zeros(1,(N/2-1)*10);
for k = 1:1:N/2-1;%5->50
    mag_Van_3(1,(bit*fsw/N)*k) = mag_Van_2(1,k);%50:0.499
    mag_Vbn_3(1,(bit*fsw/N)*k) = mag_Vbn_2(1,k);
    mag_Vcn_3(1,(bit*fsw/N)*k) = mag_Vcn_2(1,k);
end
% figure;
% semilogx(mag_Van_3,'r'),title('mag Van 3');
end_time = clock;
execution_time = end_time-start_time;
ITHDavg = (ITHD_Van+ITHD_Vbn+ITHD_Vcn)/3;

figure;
subplot(4,1,1),plot(t,s1,t,s2,t,s3);
subplot(4,1,2),plot(out(1,:));
subplot(4,1,3),plot(out(2,:));
subplot(4,1,4),plot(out(3,:));