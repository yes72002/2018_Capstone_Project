clear%1/24三相cbsvpwm就三區-跟m13接近
close all
clc
start_time=clock;
time = 0.1;
fsw = 5000;
bit = 400;
amp = 0.5;%0.2472,0.4,0.6472
fraction3=0.5;
f = 50;
n = fsw*time;
t = 0:1/fsw:time-1/fsw; 
s1 = amp*cos(2*pi*f*t);
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3);
a = 2*pi/3;
mapping=2/3*[cos(0) cos(a) cos(2*a);sin(0) sin(a) sin(2*a);1 1 1];
sindq=mapping*[s1;s2;s3];
sindqxy=[sindq(1,:);sindq(2,:)];
sindqx=sindqxy(1,:);
sindqy=sindqxy(2,:);
allsin=[s1;s2;s3];
sinmax=max(allsin);
sinmin=min(allsin);
ratio1=1;
ratio0=1;
w1=1-sinmax;
w0=-sinmin;
w=(w1*ratio1+w0*ratio0)/(ratio1+ratio0);
duty=allsin+[w;w;w];
[dutysort,dutycase] = sort(duty);
duty3rd(1,:)=dutysort(1,:);
duty2nd(1,:)=dutysort(2,:);
duty1st(1,:)=dutysort(3,:);
Duty1=duty1st-duty2nd;
Duty2=duty2nd-duty3rd;
Dutyz=1-Duty1-Duty2;
for i=1:fsw*time
%     X=[zeros(1,dutycase(3,i)) ones(1,4-dutycase(3,i))];
%     Y=[zeros(1,dutycase(2,i)) ones(1,4-dutycase(2,i))];
%     Z=[zeros(1,dutycase(1,i)) ones(1,4-dutycase(1,i))];
    v(dutycase(1,i),:)=[0 0 0 1];
    v(dutycase(2,i),:)=[0 0 1 1];
    v(dutycase(3,i),:)=[0 1 1 1];
    v1(:,i)=v(:,2);
    v2(:,i)=v(:,3);
end
A=[2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3];
for i=1:fsw*time
    vdq1(:,i)=mapping*A*v1(:,i);
    vdq2(:,i)=mapping*A*v2(:,i);
end
vdq1x=vdq1(1,:);
vdq1y=vdq1(2,:);
vdq2x=vdq2(1,:);
vdq2y=vdq2(2,:);
vrefx=sindqx;
vrefy=sindqy;
X1 = (vdq1x-vrefx).*Duty1;
X2 = (vdq2x-vrefx).*Duty2;
Xz = -vrefx.*Dutyz;
Y1 = (vdq1y-vrefy).*Duty1;
Y2 = (vdq2y-vrefy).*Duty2;
Yz = -vrefy.*Dutyz;
T1=Duty1;%無限跟有限精準度
T2=Duty2;
Tz=Dutyz;
for i=1:fsw*time
    %HDF0127
    P1 = 0;
    P2 = Xz(i)/2;
    P3 = Xz(i)/2+X1(i);
    P4 = Xz(i)/2+X1(i)+X2(i);
    P5 = Xz(i)+X1(i)+X2(i);
    HDFX1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2;
    HDFX2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDFX3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDFX4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*Tz(i)/2;
    R1 = 0;
    R2 = Yz(i)/2;
    R3 = Yz(i)/2+Y1(i);
    R4 = Yz(i)/2+Y1(i)+Y2(i);
    R5 = Yz(i)+Y1(i)+Y2(i);
    HDFY1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2;
    HDFY2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i);
    HDFY3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDFY4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*Tz(i)/2;
    %HDF0121
    P1 = 0;
    P2 = Xz(i);
    P3 = Xz(i)+X1(i)*fraction3;
    P4 = Xz(i)+X1(i)*fraction3+X2(i);
    P5 = Xz(i)+X1(i)+X2(i);
    HDF0121X1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121X2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)*fraction3;
    HDF0121X3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDF0121X4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T1(i)*(1-fraction3);   
    R1 = 0;
    R2 = Yz(i);
    R3 = Yz(i)+Y1(i)*fraction3;
    R4 = Yz(i)+Y1(i)*fraction3+Y2(i);
    R5 = Yz(i)+Y1(i)+Y2(i);
    HDF0121Y1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF0121Y2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)*fraction3;
    HDF0121Y3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDF0121Y4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T1(i)*(1-fraction3);
    %HDF7212
    P1 = 0;
    P2 = Xz(i);
    P3 = Xz(i)+X2(i)*fraction3;
    P4 = Xz(i)+X2(i)*fraction3+X1(i);
    P5 = Xz(i)+X2(i)+X1(i);
    HDF7212X1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF7212X2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T2(i)*fraction3;
    HDF7212X3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T1(i);
    HDF7212X4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Yz(i);
    R3 = Yz(i)+Y2(i)*fraction3;
    R4 = Yz(i)+Y2(i)*fraction3+Y1(i);
    R5 = Yz(i)+Y2(i)+Y1(i);
    HDF7212Y1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF7212Y2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T2(i)*fraction3;
    HDF7212Y3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T1(i);
    HDF7212Y4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*(1-fraction3);
end
HDFX = HDFX1+HDFX2+HDFX3+HDFX4;
HDFY = HDFY1+HDFY2+HDFY3+HDFY4;
HDF = HDFX+HDFY;
HDF0121X = HDF0121X1+HDF0121X2+HDF0121X3+HDF0121X4;
HDF0121Y = HDF0121Y1+HDF0121Y2+HDF0121Y3+HDF0121Y4;
HDF0121 = HDF0121X+HDF0121Y;
HDF7212X = HDF7212X1+HDF7212X2+HDF7212X3+HDF7212X4;
HDF7212Y = HDF7212Y1+HDF7212Y2+HDF7212Y3+HDF7212Y4;
HDF7212 = HDF7212X+HDF7212Y;
XXX=X1+X2+Xz;
YYY=Y1+Y2+Yz;
XXXYYY=[X1;X2;Xz;Y1;Y2;Yz];

threeHDF=[HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
allarea3=sum(abs(HDF-minHDFvalue3));

% minHDFcase3(61)=2;
% minHDFcase3(1261)=2;
% minHDFcase3(2341)=2;
% minHDFcase3(3061)=3;
% minHDFcase3(3421)=2;

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
Tzbit=bit/2-T1bit-T2bit;
TTT=T1bit+T2bit+Tzbit;
minHDFcase3=ones(1,fsw*time);
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
cbsvpwm=out;
out_2=circshift(out,[0 1]);
out3=[out;out_2];
diff1=sum(abs(out(1,:)-out_2(1,:)));
diff2=sum(abs(out(2,:)-out_2(2,:)));
diff3=sum(abs(out(3,:)-out_2(3,:)));
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
simulate=[THD_Van;ITHD_Van];
disperence=[point3;allarea3;THD_Van;ITHD_Van;ITHDavg];

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
% plot(angle*180/pi,minHDFvalue3,'m'),hold on;
legend('HDF','HDF0121','HDF7212');