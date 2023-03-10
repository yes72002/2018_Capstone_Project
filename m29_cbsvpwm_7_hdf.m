clear%1/28七相cbsvpwm舊三區-跟m21一樣
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.5;
fraction3=0.5;
f = 50;
n=fsw*time;
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
allsin=[s1;s2;s3;s4;s5;s6;s7];
sinmax=max(allsin);
sinmin=min(allsin);
ratio1=1;
ratio0=1;
w1=1-sinmax;
w0=-sinmin;
w=(w1*ratio1+w0*ratio0)/(ratio1+ratio0);
duty=allsin+[w;w;w;w;w;w;w];
[dutysort,dutycase] = sort(duty);
duty7th(1,:)=dutysort(1,:);
duty6th(1,:)=dutysort(2,:);
duty5th(1,:)=dutysort(3,:);
duty4th(1,:)=dutysort(4,:);
duty3rd(1,:)=dutysort(5,:);
duty2nd(1,:)=dutysort(6,:);
duty1st(1,:)=dutysort(7,:);
Duty1=duty1st-duty2nd;
Duty2=duty2nd-duty3rd;
Duty3=duty3rd-duty4th;
Duty4=duty4th-duty5th;
Duty5=duty5th-duty6th;
Duty6=duty6th-duty7th;
Dutyz=1-Duty1-Duty2-Duty3-Duty4-Duty5-Duty6;
for i=1:fsw*time
    v(dutycase(1,i),:)=[0 0 0 0 0 0 0 1];
    v(dutycase(2,i),:)=[0 0 0 0 0 0 1 1];
    v(dutycase(3,i),:)=[0 0 0 0 0 1 1 1];
    v(dutycase(4,i),:)=[0 0 0 0 1 1 1 1];
    v(dutycase(5,i),:)=[0 0 0 1 1 1 1 1];
    v(dutycase(6,i),:)=[0 0 1 1 1 1 1 1];
    v(dutycase(7,i),:)=[0 1 1 1 1 1 1 1];
    v1(:,i)=v(:,2);
    v2(:,i)=v(:,3);
    v3(:,i)=v(:,4);
    v4(:,i)=v(:,5);
    v5(:,i)=v(:,6);
    v6(:,i)=v(:,7);
end
A = [6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7;...
    -1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7];
for i=1:fsw*time
    vdq1(:,i)=mapping*A*v1(:,i);
    vdq2(:,i)=mapping*A*v2(:,i);
    vdq3(:,i)=mapping*A*v3(:,i);
    vdq4(:,i)=mapping*A*v4(:,i);
    vdq5(:,i)=mapping*A*v5(:,i);
    vdq6(:,i)=mapping*A*v6(:,i);
end
vdq1x=vdq1(1,:);
vdq1y=vdq1(2,:);
vdq2x=vdq2(1,:);
vdq2y=vdq2(2,:);
vdq3x=vdq3(1,:);
vdq3y=vdq3(2,:);
vdq4x=vdq4(1,:);
vdq4y=vdq4(2,:);
vdq5x=vdq5(1,:);
vdq5y=vdq5(2,:);
vdq6x=vdq6(1,:);
vdq6y=vdq6(2,:);
vrefx=sindq1x;
vrefy=sindq1y;
X1 = (vdq1x-vrefx).*Duty1;
X2 = (vdq2x-vrefx).*Duty2;
X3 = (vdq3x-vrefx).*Duty3;
X4 = (vdq4x-vrefx).*Duty4;
X5 = (vdq5x-vrefx).*Duty5;
X6 = (vdq6x-vrefx).*Duty6;
Xz = -vrefx.*Dutyz;
Y1 = (vdq1y-vrefy).*Duty1;
Y2 = (vdq2y-vrefy).*Duty2;
Y3 = (vdq3y-vrefy).*Duty3;
Y4 = (vdq4y-vrefy).*Duty4;
Y5 = (vdq5y-vrefy).*Duty5;
Y6 = (vdq6y-vrefy).*Duty6;
Yz = -vrefy.*Dutyz;
T1=Duty1;%無限跟有限精準度
T2=Duty2;
T3=Duty3;
T4=Duty4;
T5=Duty5;
T6=Duty6;
Tz=Dutyz;
for i=1:1:fsw*time%01234567
    P1=0;
    P2=Xz(i)/2;
    P3=Xz(i)/2+X1(i);
    P4=Xz(i)/2+X1(i)+X2(i);
    P5=Xz(i)/2+X1(i)+X2(i)+X3(i);
    P6=Xz(i)/2+X1(i)+X2(i)+X3(i)+X4(i);
    P7=Xz(i)/2+X1(i)+X2(i)+X3(i)+X4(i)+X5(i);
    P8=Xz(i)/2+X1(i)+X2(i)+X3(i)+X4(i)+X5(i)+X6(i);
    P9=Xz(i)+X1(i)+X2(i)+X3(i)+X4(i)+X5(i)+X6(i);
    HDFX1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2;
    HDFX2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDFX3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDFX4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T3(i);
    HDFX5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T4(i);
    HDFX6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T5(i);
    HDFX7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T6(i);
    HDFX8(i)=((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Yz(i)/2;
    R3=Yz(i)/2+Y1(i);
    R4=Yz(i)/2+Y1(i)+Y2(i);
    R5=Yz(i)/2+Y1(i)+Y2(i)+Y3(i);
    R6=Yz(i)/2+Y1(i)+Y2(i)+Y3(i)+Y4(i);
    R7=Yz(i)/2+Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i);
    R8=Yz(i)/2+Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i)+Y6(i);
    R9=Yz(i)+Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i)+Y6(i);
    HDFY1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2;
    HDFY2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i);
    HDFY3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDFY4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T3(i);
    HDFY5(i)=((R5)^2+(R5)*(R6)+(R6)^2)*T4(i);
    HDFY6(i)=((R6)^2+(R6)*(R7)+(R7)^2)*T5(i);
    HDFY7(i)=((R7)^2+(R7)*(R8)+(R8)^2)*T6(i);
    HDFY8(i)=((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Xz(i);
    P3=Xz(i)+X1(i);
    P4=Xz(i)+X1(i)+X2(i);
    P5=Xz(i)+X1(i)+X2(i)+X3(i);
    P6=Xz(i)+X1(i)+X2(i)+X3(i)+X4(i);
    P7=Xz(i)+X1(i)+X2(i)+X3(i)+X4(i)+X5(i)*fraction3;
    P8=Xz(i)+X1(i)+X2(i)+X3(i)+X4(i)+X5(i)*fraction3+X6(i);
    P9=Xz(i)+X1(i)+X2(i)+X3(i)+X4(i)+X5(i)+X6(i);
    HDF0121X1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF0121X2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T1(i);
    HDF0121X3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T2(i);
    HDF0121X4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T3(i);
    HDF0121X5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T4(i);
    HDF0121X6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3;
    HDF0121X7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T6(i);
    HDF0121X8(i)=((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Yz(i);
    R3=Yz(i)+Y1(i);
    R4=Yz(i)+Y1(i)+Y2(i);
    R5=Yz(i)+Y1(i)+Y2(i)+Y3(i);
    R6=Yz(i)+Y1(i)+Y2(i)+Y3(i)+Y4(i);
    R7=Yz(i)+Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i)*fraction3;
    R8=Yz(i)+Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i)*fraction3+Y6(i);
    R9=Yz(i)+Y1(i)+Y2(i)+Y3(i)+Y4(i)+Y5(i)+Y6(i);
    HDF0121Y1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF0121Y2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T1(i);
    HDF0121Y3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T2(i);
    HDF0121Y4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T3(i);
    HDF0121Y5(i)=((R5)^2+(R5)*(R6)+(R6)^2)*T4(i);
    HDF0121Y6(i)=((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3;
    HDF0121Y7(i)=((R7)^2+(R7)*(R8)+(R8)^2)*T6(i);
    HDF0121Y8(i)=((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Xz(i);
    P3=Xz(i)+X6(i);
    P4=Xz(i)+X6(i)+X5(i);
    P5=Xz(i)+X6(i)+X5(i)+X4(i);
    P6=Xz(i)+X6(i)+X5(i)+X4(i)+X3(i);
    P7=Xz(i)+X6(i)+X5(i)+X4(i)+X3(i)+X2(i)*fraction3;
    P8=Xz(i)+X6(i)+X5(i)+X4(i)+X3(i)+X2(i)*fraction3+X1(i);
    P9=Xz(i)+X6(i)+X5(i)+X4(i)+X3(i)+X2(i)+X1(i);
    HDF7212X1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i);
    HDF7212X2(i)=((P2)^2+(P2)*(P3)+(P3)^2)*T6(i);
    HDF7212X3(i)=((P3)^2+(P3)*(P4)+(P4)^2)*T5(i);
    HDF7212X4(i)=((P4)^2+(P4)*(P5)+(P5)^2)*T4(i);
    HDF7212X5(i)=((P5)^2+(P5)*(P6)+(P6)^2)*T3(i);
    HDF7212X6(i)=((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3;
    HDF7212X7(i)=((P7)^2+(P7)*(P8)+(P8)^2)*T1(i);
    HDF7212X8(i)=((P8)^2+(P8)*(P1)+(P1)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Yz(i);
    R3=Yz(i)+Y6(i);
    R4=Yz(i)+Y6(i)+Y5(i);
    R5=Yz(i)+Y6(i)+Y5(i)+Y4(i);
    R6=Yz(i)+Y6(i)+Y5(i)+Y4(i)+Y3(i);
    R7=Yz(i)+Y6(i)+Y5(i)+Y4(i)+Y3(i)+Y2(i)*fraction3;
    R8=Yz(i)+Y6(i)+Y5(i)+Y4(i)+Y3(i)+Y2(i)*fraction3+Y1(i);
    R9=Yz(i)+Y6(i)+Y5(i)+Y4(i)+Y3(i)+Y2(i)+Y1(i);
    HDF7212Y1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i);
    HDF7212Y2(i)=((R2)^2+(R2)*(R3)+(R3)^2)*T6(i);
    HDF7212Y3(i)=((R3)^2+(R3)*(R4)+(R4)^2)*T5(i);
    HDF7212Y4(i)=((R4)^2+(R4)*(R5)+(R5)^2)*T4(i);
    HDF7212Y5(i)=((R5)^2+(R5)*(R6)+(R6)^2)*T3(i);
    HDF7212Y6(i)=((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3;
    HDF7212Y7(i)=((R7)^2+(R7)*(R8)+(R8)^2)*T1(i);
    HDF7212Y8(i)=((R8)^2+(R8)*(R1)+(R1)^2)*T2(i)*(1-fraction3);
end
HDFX=HDFX1+HDFX2+HDFX3+HDFX4+HDFX5+HDFX6+HDFX7+HDFX8;
HDFY=HDFY1+HDFY2+HDFY3+HDFY4+HDFY5+HDFY6+HDFY7+HDFY8;
HDF=HDFX+HDFY;
HDF0121X=HDF0121X1+HDF0121X2+HDF0121X3+HDF0121X4+HDF0121X5+HDF0121X6+HDF0121X7+HDF0121X8;
HDF0121Y=HDF0121Y1+HDF0121Y2+HDF0121Y3+HDF0121Y4+HDF0121Y5+HDF0121Y6+HDF0121Y7+HDF0121Y8;
HDF0121=HDF0121X+HDF0121Y;
HDF7212X=HDF7212X1+HDF7212X2+HDF7212X3+HDF7212X4+HDF7212X5+HDF7212X6+HDF7212X7+HDF7212X8;
HDF7212Y=HDF7212Y1+HDF7212Y2+HDF7212Y3+HDF7212Y4+HDF7212Y5+HDF7212Y6+HDF7212Y7+HDF7212Y8;
HDF7212=HDF7212X+HDF7212Y;

threeHDF=[HDF;HDF0121;HDF7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
allarea3=sum(abs(HDF-minHDFvalue3));

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
T5bit=round(T5*bit/2);
T6bit=round(T6*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit-T5bit-T6bit;
for i=1:fsw*time%主要的大程式碼
   switch minHDFcase3(i)
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
fprintf('out1切換次數=%f\n',diff1);
fprintf('out2切換次數=%f\n',diff2);
fprintf('out3切換次數=%f\n',diff3);
fprintf('out4切換次數=%f\n',diff4);
fprintf('out5切換次數=%f\n',diff5);
fprintf('out6切換次數=%f\n',diff6);
fprintf('out7切換次數=%f\n',diff7);

N = time*fsw*bit;
Van=[6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vbn=[-1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vcn=[-1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7]*out;
Vdn=[-1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7]*out;
Ven=[-1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7]*out;
Vfn=[-1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7]*out;
Vgn=[-1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7]*out;
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
fft_Vfn=fft(Vfn,N)/N;
mag_Vfn=abs(fft_Vfn)*2;
fft_Vgn=fft(Vgn,N)/N;
mag_Vgn=abs(fft_Vgn)*2;
k=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:到中間
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
mag_Vfn_1=mag_Vfn(1,k);
mag_Vgn_1=mag_Vgn(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(倍頻^2/基頻^2)
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
%ITHD=頻率除以本身值
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
mag_ithd_out6=zeros(1,N/5/2-1);
mag_ithd_out7=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%末減首除以公差+1=[(N/2-4)-6]/5+1
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
%該算的,此行以上都算出來了,接下來是畫圖
n2=2:1:N/2;%去頭(頭第一個值有誤為零)
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
disperence=[point3;allarea3;THD_Van;ITHD_Van;ITHD_Vbn;ITHDavg];
figure;%allHDF
angle=atan2(sindq1y,sindq1x);
for i=1:1:fsw*time
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
plot(angle*180/pi,HDF0121,'b'),hold on;
plot(angle*180/pi,HDF7212,'g'),hold on;
plot(angle*180/pi,minHDFvalue3,'m'),hold on;