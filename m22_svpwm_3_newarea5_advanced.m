clear%2/6穝き跋纔て
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
s2 = amp*cos(2*pi*f*t-2*pi/3); 
s3 = amp*cos(2*pi*f*t-4*pi/3); 
a = 2*pi/3;
mapping=2/3*[cos(0) cos(a) cos(2*a);sin(0) sin(a) sin(2*a);1 1 1];
sindqxy=mapping*[s1;s2;s3];
sindqx=sindqxy(1,:);
sindqy=sindqxy(2,:);
angle=atan2(sindqy,sindqx);
for i=1:1:fsw*time
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
A = [2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3];
for i=1:6
    binary=dec2bin(i,3);
    bin=binary-48;
    DQ(:,i)=mapping*A*[bin(1);bin(2);bin(3)];
    M1(:,i)=[DQ(1,i);DQ(2,i)];
end
% figure,plot(sindq1y,sindq1x),hold on;
% quiver(zeros(1,126),zeros(1,126),DQ(1,:),DQ(2,:)),title('dq-plane 1'); 
% for i=1:126
%     k=num2str(i);
%     text(DQ(1,i),DQ(2,i),k);
% end
for i=1:fsw*time%耞sector
    if angle(i)>0 && angle(i)<=pi/3;%sector1
        sector(i)=1;
        sb1(i)=4;
        sb2(i)=6;
    end
    if angle(i)>pi/3 && angle(i)<=2*pi/3;%sector2
        sector(i)=2;
        sb1(i)=2;
        sb2(i)=6;
    end   
    if angle(i)>2*pi/3 && angle(i)<=3*pi/3;%sector3
        sector(i)=3;
        sb1(i)=2;
        sb2(i)=3;
    end
    if angle(i)>3*pi/3 && angle(i)<=4*pi/3;%sector4
        sector(i)=4;
        sb1(i)=1;
        sb2(i)=3;
    end
    if angle(i)>4*pi/3 && angle(i)<=5*pi/3;%sector5
        sector(i)=5;
        sb1(i)=1;
        sb2(i)=5;
    end 
    if angle(i)>5*pi/3 && angle(i)<=6*pi/3;%sector6
        sector(i)=6;
        sb1(i)=4;
        sb2(i)=5;
    end
    duty(:,i)=[M1(:,sb1(i)) M1(:,sb2(i))]\[sindqx(i);sindqy(i)];
    v1(:,i)=dec2bin(sb1(i),3)-48;
    v2(:,i)=dec2bin(sb2(i),3)-48;
end
T1=duty(1,:);%礚&Τ弘非
T2=duty(2,:);
Tz=1-T1-T2;
for i=1:1:fsw*time%玻ネQ1,Q2,Qz
    if mod(sector(i),2)==1%1,sector碞琌计0,sector碞琌案计
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
QQQ_p1=Qz+Q1+Q2;
DDD_p2=Dz+D1+D2;
for i=1:fsw*time%Plane1
    P1 = 0;%0127
    P2 = Qz(i)/2;
    P3 = Qz(i)/2+Q1(i);
    P4 = Qz(i)/2+Q1(i)+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i);
    HDFQ(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
           +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
           +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
           +((P4)^2+(P4)*(P5)+(P5)^2)*Tz(i)/2;
    R1 = 0;
    R2 = Dz(i)/2;
    R3 = Dz(i)/2+D1(i);
    R4 = Dz(i)/2+D1(i)+D2(i);
    R5 = Dz(i)+D1(i)+D2(i);
    HDFD(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
           +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
           +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
           +((R4)^2+(R4)*(R5)+(R5)^2)*Tz(i)/2;
    P1 = 0;%0121
    P2 = Qz(i);
    P3 = Qz(i)+Q1(i)*fraction3;
    P4 = Qz(i)+Q1(i)*fraction3+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i);
    HDF0121Q(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)*fraction3...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T1(i)*(1-fraction3);    
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D1(i)*fraction3;
    R4 = Dz(i)+D1(i)*fraction3+D2(i);
    R5 = Dz(i)+D1(i)+D2(i);
    HDF0121D(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)*fraction3...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T1(i)*(1-fraction3);
    P1 = 0;%7212
    P2 = Qz(i);
    P3 = Qz(i)+Q2(i)*fraction3;
    P4 = Qz(i)+Q2(i)*fraction3+Q1(i);
    P5 = Qz(i)+Q2(i)+Q1(i);
    HDF7212Q(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T2(i)*fraction3...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T1(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D2(i)*fraction3;
    R4 = Dz(i)+D2(i)*fraction3+D1(i);
    R5 = Dz(i)+D2(i)+D1(i);
    HDF7212D(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T2(i)*fraction3...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T1(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*(1-fraction3);
    P1 = 0;%new0121
    P2 = Qz(i);
    P3 = Qz(i)+Q1(i)*fraction5;
    P4 = Qz(i)+Q1(i)*fraction5+Q2(i);
    P5 = Qz(i)+Q1(i)+Q2(i);
    HDFnew0121Q(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)*fraction5...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T1(i)*(1-fraction5);    
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D1(i)*fraction5;
    R4 = Dz(i)+D1(i)*fraction5+D2(i);
    R5 = Dz(i)+D1(i)+D2(i);
    HDFnew0121D(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)*fraction5...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T1(i)*(1-fraction5);
    P1 = 0;%new7212
    P2 = Qz(i);
    P3 = Qz(i)+Q2(i)*fraction5;
    P4 = Qz(i)+Q2(i)*fraction5+Q1(i);
    P5 = Qz(i)+Q2(i)+Q1(i);
    HDFnew7212Q(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T2(i)*fraction5...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T1(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*(1-fraction5);
    R1 = 0;
    R2 = Dz(i);
    R3 = Dz(i)+D2(i)*fraction5;
    R4 = Dz(i)+D2(i)*fraction5+D1(i);
    R5 = Dz(i)+D2(i)+D1(i);
    HDFnew7212D(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T2(i)*fraction5...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T1(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*(1-fraction5);
end
HDF = HDFQ + HDFD;
HDF0121 = HDF0121Q + HDF0121D;
HDF7212 = HDF7212Q + HDF7212D;
HDFnew0121 = HDFnew0121Q + HDFnew0121D;
HDFnew7212 = HDFnew7212Q + HDFnew7212D;
threeHDF = [HDF;HDF0121;HDF7212];
fiveHDF = [HDF;HDF0121;HDF7212;HDFnew0121;HDFnew7212];
[minHDFvalue3,minHDFcase3] = min(threeHDF);
[minHDFvalue5,minHDFcase5] = min(fiveHDF);
point3 = sum(minHDFcase3==2,2) + sum(minHDFcase3==3,2);
point5 = sum(minHDFcase5==4,2) + sum(minHDFcase5==5,2);
allarea3 = sum(abs(HDF-minHDFvalue3));
allarea5 = sum(abs(minHDFvalue3-minHDFvalue5));
T1bit = round(T1*bit/2);
T2bit = round(T2*bit/2);
Tzbit  =bit/2 - T1bit - T2bit;
for i=1:fsw*time%璶祘Α絏
   switch minHDFcase5(i)
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
       case 4%new0121
           out000=repmat([0;0;0],1,Tzbit(i));
           out001=repmat(v1(:,i),1,floor(T1bit(i)*fraction5));
           out011=repmat(v2(:,i),1,T2bit(i));
           out111=repmat(v1(:,i),1,T1bit(i)-floor(T1bit(i)*fraction5));
           outnew0121=[out000 out001 out011 out111 out111 out011 out001 out000];
           out(:,bit*(i-1)+1:bit*i)=outnew0121;
       case 5%new7212
           out000=repmat([1;1;1],1,Tzbit(i));
           out001=repmat(v2(:,i),1,floor(T2bit(i)*fraction5));
           out011=repmat(v1(:,i),1,T1bit(i));
           out111=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)*fraction5));
           outnew7212=[out000 out001 out011 out111 out111 out011 out001 out000];
           out(:,bit*(i-1)+1:bit*i)=outnew7212;
   end
end
out_2=circshift(out,[0 1]);
out_3=[out;out_2];
diff1=sum(abs(out(1,:)-out_2(1,:)));
diff2=sum(abs(out(2,:)-out_2(2,:)));
diff3=sum(abs(out(3,:)-out_2(3,:)));
fprintf('out1ち传Ω计=%f\n',diff1);
fprintf('out2ち传Ω计=%f\n',diff2);
fprintf('out3ち传Ω计=%f\n',diff3);

N = time*fsw*bit;
Van=[ 2/3 -1/3 -1/3]*out;
Vbn=[-1/3  2/3 -1/3]*out;
Vcn=[-1/3 -1/3  2/3]*out;
fft_Van=fft(Van,N)/N;%6:0.499ㄢ狠常Τ
mag_Van=abs(fft_Van)*2;
fft_Vbn=fft(Vbn,N)/N;
mag_Vbn=abs(fft_Vbn)*2;
fft_Vcn=fft(Vcn,N)/N;
mag_Vcn=abs(fft_Vcn)*2;
k=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:い丁
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
THD_Van=100*((sum(mag_Van_1.^2)-(mag_Van_1(1,1).^2))/(mag_Van_1(1,1).^2)).^(1/2);
THD_Vbn=100*((sum(mag_Vbn_1.^2)-(mag_Vbn_1(1,1).^2))/(mag_Vbn_1(1,1).^2)).^(1/2);
THD_Vcn=100*((sum(mag_Vcn_1.^2)-(mag_Vcn_1(1,1).^2))/(mag_Vcn_1(1,1).^2)).^(1/2);
fprintf('VTHD1=%f\n',THD_Van);
fprintf('VTHD2=%f\n',THD_Vbn);
fprintf('VTHD3=%f\n',THD_Vcn);
for i=1:length(mag_Van_1)
    mag_ithd_out1(i)=mag_Van_1(i)/(10*i); 
    mag_ithd_out2(i)=mag_Vbn_1(i)/(10*i); 
    mag_ithd_out3(i)=mag_Vcn_1(i)/(10*i);  
end
ITHD_Van=100*((sum(mag_ithd_out1.^2)-(mag_ithd_out1(1,1).^2))/(mag_ithd_out1(1,1).^2)).^(1/2);
ITHD_Vbn=100*((sum(mag_ithd_out2.^2)-(mag_ithd_out2(1,1).^2))/(mag_ithd_out2(1,1).^2)).^(1/2);
ITHD_Vcn=100*((sum(mag_ithd_out3.^2)-(mag_ithd_out3(1,1).^2))/(mag_ithd_out3(1,1).^2)).^(1/2);
fprintf('ITHD1=%f\n',ITHD_Van);
fprintf('ITHD2=%f\n',ITHD_Vbn);
fprintf('ITHD3=%f\n',ITHD_Vcn);
end_time=clock;
execution_time=end_time-start_time;
ITHDavg=(ITHD_Van+ITHD_Vbn+ITHD_Vcn)/3;
disperence=[point5;allarea5;THD_Van;ITHD_Van;ITHD_Vbn;ITHDavg];
% figure;%allHDF
% plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
% plot(angle*180/pi,HDF0121,'b'),hold on;
% plot(angle*180/pi,HDF7212,'g'),hold on;
% plot(angle*180/pi,minHDFvalue3,'m'),hold on;
% plot(angle*180/pi,HDFnew0121,'c'),hold on;
% plot(angle*180/pi,HDFnew7212,'y'),hold on;
% plot(angle*180/pi,minHDFvalue5,'k'),hold on;