clear%4/18三相
close all
clc
start_time = clock;
time  =  0.1;
fsw  =  36000;
bit  =  400;
amp  =  0.4;
fraction3 = 0.5;
f  =  50;
n  =  fsw*time;
t  =  0:1/fsw:time-1/fsw;
s1  =  amp*cos(2*pi*f*t);
s2  =  amp*cos(2*pi*f*t-2*pi/3);
s3  =  amp*cos(2*pi*f*t-4*pi/3);
a = 2*pi/3;
mapping = 2/3*[cos(0) cos(a) cos(2*a);sin(0) sin(a) sin(2*a);1 1 1];
sindq = mapping*[s1;s2;s3];
sindqxy = [sindq(1,:);sindq(2,:)];
sindqx = sindqxy(1,:);
sindqy = sindqxy(2,:);
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
updown = floor(allsin);
for i = 1:fsw*time
    v(dutycase(1,i),:) = [0 0 0 1];
    v(dutycase(2,i),:) = [0 0 1 1];
    v(dutycase(3,i),:) = [0 1 1 1];
    v0(:,i) = v(:,1);
    v1(:,i) = v(:,2);
    v2(:,i) = v(:,3);
    v3(:,i) = v(:,4);
    ud = updown(:,i);
    v0real(:,i) = v0(:,i)+ud;
    v1real(:,i) = v1(:,i)+ud;
    v2real(:,i) = v2(:,i)+ud;
    v3real(:,i) = v3(:,i)+ud;
    A = [2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3];
    vdq0(:,i) = mapping*A*v0real(:,i);
    vdq1(:,i) = mapping*A*v1real(:,i);
    vdq2(:,i) = mapping*A*v2real(:,i);
    vdq3(:,i) = mapping*A*v3real(:,i);
end
vdq0x=vdq0(1,:);
vdq0y=vdq0(2,:);
vdq1x=vdq1(1,:);
vdq1y=vdq1(2,:);
vdq2x=vdq2(1,:);
vdq2y=vdq2(2,:);
vdq3x=vdq3(1,:);
vdq3y=vdq3(2,:);
vrefx=sindqx;
vrefy=sindqy;
X0 = (vdq0x-vrefx).*Dutyz/2;
X1 = (vdq1x-vrefx).*Duty1;
X2 = (vdq2x-vrefx).*Duty2;
X3 = (vdq3x-vrefx).*Dutyz/2;
Xz = X0 + X3;
Y0 = (vdq0y-vrefy).*Dutyz/2;
Y1 = (vdq1y-vrefy).*Duty1;
Y2 = (vdq2y-vrefy).*Duty2;
Y3 = (vdq3y-vrefy).*Dutyz/2;
Yz = Y0 + Y3;

T1 = Duty1;%無限跟有限精準度
T2 = Duty2;
Tz = Dutyz;
for i=1:fsw*time
    P1 = 0;%HDF0127
    P2 = X0(i);
    P3 = X0(i)+X1(i);
    P4 = X0(i)+X1(i)+X2(i);
    P5 = X0(i)+X1(i)+X2(i)+X3(i);
    HDFX(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
          +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
          +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
          +((P4)^2+(P4)*(P5)+(P5)^2)*Tz(i)/2;
    R1 = 0;
    R2 = Y0(i);
    R3 = Y0(i)+Y1(i);
    R4 = Y0(i)+Y1(i)+Y2(i);
    R5 = Y0(i)+Y1(i)+Y2(i)+Y3(i);
    HDFY(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
          +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
          +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
          +((R4)^2+(R4)*(R5)+(R5)^2)*Tz(i)/2;
    P1 = 0;%HDF0121
    P2 = Xz(i);
    P3 = Xz(i)+X1(i)*fraction3;
    P4 = Xz(i)+X1(i)*fraction3+X2(i);
    P5 = Xz(i)+X1(i)+X2(i);
    HDF0121X(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
             +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)*fraction3...
             +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
             +((P4)^2+(P4)*(P5)+(P5)^2)*T1(i)*(1-fraction3);   
    R1 = 0;
    R2 = Yz(i);
    R3 = Yz(i)+Y1(i)*fraction3;
    R4 = Yz(i)+Y1(i)*fraction3+Y2(i);
    R5 = Yz(i)+Y1(i)+Y2(i);
    HDF0121Y(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
             +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)*fraction3...
             +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
             +((R4)^2+(R4)*(R5)+(R5)^2)*T1(i)*(1-fraction3);
    P1 = 0;%HDF7212
    P2 = Xz(i);
    P3 = Xz(i)+X2(i)*fraction3;
    P4 = Xz(i)+X2(i)*fraction3+X1(i);
    P5 = Xz(i)+X2(i)+X1(i);
    HDF7212X(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
             +((P2)^2+(P2)*(P3)+(P3)^2)*T2(i)*fraction3...
             +((P3)^2+(P3)*(P4)+(P4)^2)*T1(i)...
             +((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Yz(i);
    R3 = Yz(i)+Y2(i)*fraction3;
    R4 = Yz(i)+Y2(i)*fraction3+Y1(i);
    R5 = Yz(i)+Y2(i)+Y1(i);
    HDF7212Y(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
             +((R2)^2+(R2)*(R3)+(R3)^2)*T2(i)*fraction3...
             +((R3)^2+(R3)*(R4)+(R4)^2)*T1(i)...
             +((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*(1-fraction3);
end
HDF = HDFX + HDFY;
HDF0121 = HDF0121X + HDF0121Y;
HDF7212 = HDF7212X + HDF7212Y;
Fdist=sqrt(sum(HDF));
Fdist0121=sqrt(sum(HDF0121));
Fdist7212=sqrt(sum(HDF7212));
fprintf('Fdist = %f\n',Fdist);
fprintf('Fdist0121 = %f\n',Fdist0121);
fprintf('Fdist7212 = %f\n',Fdist7212);

threeHDF = [HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3] = min(threeHDF);
point3 = sum(minHDFcase3==2,2) + sum(minHDFcase3==3,2);
allarea3 = sum(abs(HDF-minHDFvalue3));

T1bit = round(T1*bit/2);
T2bit = round(T2*bit/2);
Tzbit = bit/2-T1bit-T2bit;
updown = floor(allsin);
minHDFcase3=ones(1,fsw*time)+2;
for i = 1:fsw*time%主要的大程式碼
    ud=updown(:,i);
    switch minHDFcase3(i)
        case 1%0127
            out000 = repmat([0;0;0]+ud,1,floor(Tzbit(i)/2));
            out001 = repmat(v1(:,i)+ud,1,T1bit(i));
            out011 = repmat(v2(:,i)+ud,1,T2bit(i));
            out111 = repmat([1;1;1]+ud,1,Tzbit(i)-floor(Tzbit(i)/2));
            out0127 = [out000 out001 out011 out111 out111 out011 out001 out000];
            out(:,bit*(i-1)+1:bit*i) = out0127;
        case 2%0121
            out000 = repmat([0;0;0]+ud,1,Tzbit(i));
            out001 = repmat(v1(:,i)+ud,1,floor(T1bit(i)/2));
            out011 = repmat(v2(:,i)+ud,1,T2bit(i));
            out111 = repmat(v1(:,i)+ud,1,T1bit(i)-floor(T1bit(i)/2));
            out0121 = [out000 out001 out011 out111 out111 out011 out001 out000];
            out(:,bit*(i-1)+1:bit*i) = out0121;
        case 3%7212
            out000 = repmat([1;1;1]+ud,1,Tzbit(i));
            out001 = repmat(v2(:,i)+ud,1,floor(T2bit(i)/2));
            out011 = repmat(v1(:,i)+ud,1,T1bit(i));
            out111 = repmat(v2(:,i)+ud,1,T2bit(i)-floor(T2bit(i)/2));
            out7212 = [out000 out001 out011 out111 out111 out011 out001 out000];
            out(:,bit*(i-1)+1:bit*i) = out7212;
   end
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
simulate = [THD_Van;ITHD_Van];

figure;
subplot(4,1,1),plot(t,s1,t,s2,t,s3);
subplot(4,1,2),plot(out(1,:));
subplot(4,1,3),plot(out(2,:));
subplot(4,1,4),plot(out(3,:));

figure;%allHDF
angle=atan2(sindqy,sindqx);
for i=1:1:n
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
plot(angle*180/pi,HDF0121,'b'),hold on;
plot(angle*180/pi,HDF7212,'g'),hold on;
legend('HDF', 'HDF0121','HDF7212')