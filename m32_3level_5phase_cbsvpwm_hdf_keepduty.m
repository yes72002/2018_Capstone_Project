clear%4/29����
close all
clc
start_time = clock;
time = 0.1;
fsw = 3e4;
bit = 400;
amp = 0.9;%0.2472,0.4,0.6472
fraction3 = 0.5;
f = 50;
n = fsw*time;
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
sindq_p1x=sindqxy(1,:);
sindq_p1y=sindqxy(2,:);
sindq_p2x=sindqxy(3,:);
sindq_p2y=sindqxy(4,:);
allsin=[s1;s2;s3;s4;s5];
sinmax = max(allsin);
sinmin = min(allsin);

w=(-sinmin-sinmax)/2;
allsin=allsin+[w;w;w;w;w];

duty = allsin-floor(allsin);
[dutysort,dutycase] = sort(duty);
duty5th(1,:) = dutysort(1,:);
duty4th(1,:) = dutysort(2,:);
duty3rd(1,:) = dutysort(3,:);
duty2nd(1,:) = dutysort(4,:);
duty1st(1,:) = dutysort(5,:);
Duty0 = 1-duty1st;
Duty1 = duty1st-duty2nd;
Duty2 = duty2nd-duty3rd;
Duty3 = duty3rd-duty4th;
Duty4 = duty4th-duty5th;
Duty5 = duty5th;
updown = floor(allsin);
for i = 1:fsw*time
    v(dutycase(1,i),:) = [0 0 0 0 0 1];
    v(dutycase(2,i),:) = [0 0 0 0 1 1];
    v(dutycase(3,i),:) = [0 0 0 1 1 1];
    v(dutycase(4,i),:) = [0 0 1 1 1 1];
    v(dutycase(5,i),:) = [0 1 1 1 1 1];
    v0(:,i) = v(:,1);
    v1(:,i) = v(:,2);
    v2(:,i) = v(:,3);
    v3(:,i) = v(:,4);
    v4(:,i) = v(:,5);
    v5(:,i) = v(:,6);
    ud = updown(:,i);
    A = [4/5 -1/5 -1/5 -1/5 -1/5;...
    -1/5 4/5 -1/5 -1/5 -1/5;...
    -1/5 -1/5 4/5 -1/5 -1/5;...
    -1/5 -1/5 -1/5 4/5 -1/5;...
    -1/5 -1/5 -1/5 -1/5 4/5];
    vdq0(:,i) = mapping*A*(v0(:,i) + ud);
    vdq1(:,i) = mapping*A*(v1(:,i) + ud);
    vdq2(:,i) = mapping*A*(v2(:,i) + ud);
    vdq3(:,i) = mapping*A*(v3(:,i) + ud);
    vdq4(:,i) = mapping*A*(v4(:,i) + ud);
    vdq5(:,i) = mapping*A*(v5(:,i) + ud);
end
X0_p1 = (vdq0(1,:)-sindq_p1x).*Duty0;%v0�bdq����1�Wx�b
X1_p1 = (vdq1(1,:)-sindq_p1x).*Duty1;%v1�bdq����1�Wx�b
X2_p1 = (vdq2(1,:)-sindq_p1x).*Duty2;%v2�bdq����1�Wx�b
X3_p1 = (vdq3(1,:)-sindq_p1x).*Duty3;%v3�bdq����1�Wx�b
X4_p1 = (vdq4(1,:)-sindq_p1x).*Duty4;%v4�bdq����1�Wx�b
X5_p1 = (vdq5(1,:)-sindq_p1x).*Duty5;%v5�bdq����1�Wx�b
Y0_p1 = (vdq0(2,:)-sindq_p1y).*Duty0;%v1�bdq����1�Wy�b
Y1_p1 = (vdq1(2,:)-sindq_p1y).*Duty1;%v1�bdq����1�Wy�b
Y2_p1 = (vdq2(2,:)-sindq_p1y).*Duty2;%v2�bdq����1�Wy�b
Y3_p1 = (vdq3(2,:)-sindq_p1y).*Duty3;%v3�bdq����1�Wy�b
Y4_p1 = (vdq4(2,:)-sindq_p1y).*Duty4;%v4�bdq����1�Wy�b
Y5_p1 = (vdq5(2,:)-sindq_p1y).*Duty5;%v5�bdq����1�Wy�b
X0_p2 = (vdq0(3,:)-sindq_p2x).*Duty0;%v1�bdq����2�Wx�b
X1_p2 = (vdq1(3,:)-sindq_p2x).*Duty1;%v1�bdq����2�Wx�b
X2_p2 = (vdq2(3,:)-sindq_p2x).*Duty2;%v2�bdq����2�Wx�b
X3_p2 = (vdq3(3,:)-sindq_p2x).*Duty3;%v3�bdq����2�Wx�b
X4_p2 = (vdq4(3,:)-sindq_p2x).*Duty4;%v4�bdq����2�Wx�b
X5_p2 = (vdq5(3,:)-sindq_p2x).*Duty5;%v5�bdq����2�Wx�b
Y0_p2 = (vdq0(4,:)-sindq_p2y).*Duty0;%v0�bdq����2�Wy�b
Y1_p2 = (vdq1(4,:)-sindq_p2y).*Duty1;%v1�bdq����2�Wy�b
Y2_p2 = (vdq2(4,:)-sindq_p2y).*Duty2;%v2�bdq����2�Wy�b
Y3_p2 = (vdq3(4,:)-sindq_p2y).*Duty3;%v3�bdq����2�Wy�b
Y4_p2 = (vdq4(4,:)-sindq_p2y).*Duty4;%v4�bdq����2�Wy�b
Y5_p2 = (vdq5(4,:)-sindq_p2y).*Duty5;%v5�bdq����2�Wy�b
Xz_p1 = X0_p1 + X5_p1;
Xz_p2 = X0_p2 + X5_p2;
Yz_p1 = Y0_p1 + Y5_p1;
Yz_p2 = Y0_p2 + Y5_p2;
T0 = Duty0;
T1 = Duty1;%�L���򦳭���ǫ�
T2 = Duty2;
T3 = Duty3;
T4 = Duty4;
T5 = Duty5;
Tz = Duty0 + Duty5;
for i=1:fsw*time%Plane1
    P1 = 0;%HDF012347
    P2 = X0_p1(i);
    P3 = X0_p1(i)+X1_p1(i);
    P4 = X0_p1(i)+X1_p1(i)+X2_p1(i);
    P5 = X0_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i);
    P6 = X0_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i);
    P7 = X0_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i);
    HDFX_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*T0(i)...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i);
    R1 = 0;
    R2 = Y0_p1(i);
    R3 = Y0_p1(i)+Y1_p1(i);
    R4 = Y0_p1(i)+Y1_p1(i)+Y2_p1(i);
    R5 = Y0_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i);
    R6 = Y0_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i);
    R7 = Y0_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i);
    HDFY_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*T0(i)...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i);
    P1 = 0;%HDF012343
    P2 = Xz_p1(i);
    P3 = Xz_p1(i)+X1_p1(i);
    P4 = Xz_p1(i)+X1_p1(i)+X2_p1(i);
    P5 = Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)*fraction3;
    P6 = Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)*fraction3+X4_p1(i);
    P7 = Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i);
    HDF0121X_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p1(i);
    R3 = Yz_p1(i)+Y1_p1(i);
    R4 = Yz_p1(i)+Y1_p1(i)+Y2_p1(i);
    R5 = Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)*fraction3;
    R6 = Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)*fraction3+Y4_p1(i);
    R7 = Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i);
    HDF0121Y_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Xz_p1(i);
    P3 = Xz_p1(i)+X4_p1(i);
    P4 = Xz_p1(i)+X4_p1(i)+X3_p1(i);
    P5 = Xz_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)*fraction3;
    P6 = Xz_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)*fraction3+X1_p1(i);
    P7 = Xz_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)+X1_p1(i);
    HDF7212X_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p1(i);
    R3 = Yz_p1(i)+Y4_p1(i);
    R4 = Yz_p1(i)+Y4_p1(i)+Y3_p1(i);
    R5 = Yz_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)*fraction3;
    R6 = Yz_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)*fraction3+Y1_p1(i);
    R7 = Yz_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)+Y1_p1(i);
    HDF7212Y_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*(1-fraction3);
end
for i=1:fsw*time%Plane2
    P1 = 0;%HDF012347
    P2 = X0_p2(i);
    P3 = X0_p2(i)+X1_p2(i);
    P4 = X0_p2(i)+X1_p2(i)+X2_p2(i);
    P5 = X0_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i);
    P6 = X0_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i);
    P7 = X0_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i);
    HDFX_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*T0(i) ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i);
    R1 = 0;
    R2 = Y0_p2(i);
    R3 = Y0_p2(i)+Y1_p2(i);
    R4 = Y0_p2(i)+Y1_p2(i)+Y2_p2(i);
    R5 = Y0_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i);
    R6 = Y0_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i);
    R7 = Y0_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i);
    HDFY_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*T0(i) ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i);
    P1 = 0;%HDF012343
    P2 = Xz_p2(i);
    P3 = Xz_p2(i)+X1_p2(i);
    P4 = Xz_p2(i)+X1_p2(i)+X2_p2(i);
    P5 = Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)*fraction3;
    P6 = Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)*fraction3+X4_p2(i);
    P7 = Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i);
    HDF0121X_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p2(i);
    R3 = Yz_p2(i)+Y1_p2(i);
    R4 = Yz_p2(i)+Y1_p2(i)+Y2_p2(i);
    R5 = Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)*fraction3;
    R6 = Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)*fraction3+Y4_p2(i);
    R7 = Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i);
    HDF0121Y_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Xz_p2(i);
    P3 = Xz_p2(i)+X4_p2(i);
    P4 = Xz_p2(i)+X4_p2(i)+X3_p2(i);
    P5 = Xz_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)*fraction3;
    P6 = Xz_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)*fraction3+X1_p2(i);
    P7 = Xz_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)+X1_p2(i);
    HDF7212X_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p2(i);
    R3 = Yz_p2(i)+Y4_p2(i);
    R4 = Yz_p2(i)+Y4_p2(i)+Y3_p2(i);
    R5 = Yz_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)*fraction3;
    R6 = Yz_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)*fraction3+Y1_p2(i);
    R7 = Yz_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)+Y1_p2(i);
    HDF7212Y_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*(1-fraction3);
end
HDF = HDFX_p1 + HDFY_p1 + HDFX_p2 + HDFY_p2;
HDF0121 = HDF0121X_p1 + HDF0121Y_p1 + HDF0121X_p2 + HDF0121Y_p2;
HDF7212 = HDF7212X_p1 + HDF7212Y_p1 + HDF7212X_p2 + HDF7212Y_p2;
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

T0bit = floor(T0*bit/2);
T1bit = round(T1*bit/2);
T2bit = round(T2*bit/2);
T3bit = round(T3*bit/2);
T4bit = round(T4*bit/2);
T5bit = bit/2-T0bit-T1bit-T2bit-T3bit-T4bit;
Tzbit = T0bit+T5bit;
updown = floor(allsin);
for i = 1:fsw*time%�D�n���j�{���X
    ud=updown(:,i);
    switch minHDFcase3(i)
        case 1%012347
            out0=repmat([0;0;0;0;0]+ud,1,T0bit(i));
            out1=repmat(v1(:,i)+ud,1,T1bit(i));
            out2=repmat(v2(:,i)+ud,1,T2bit(i));
            out3=repmat(v3(:,i)+ud,1,T3bit(i));
            out4=repmat(v4(:,i)+ud,1,T4bit(i));
            out7=repmat([1;1;1;1;1]+ud,1,T5bit(i));
            out012347=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
            out(:,bit*(i-1)+1:bit*i)=out012347;
        case 2%012343
            out0=repmat([0;0;0;0;0]+ud,1,Tzbit(i));
            out1=repmat(v1(:,i)+ud,1,T1bit(i));
            out2=repmat(v2(:,i)+ud,1,T2bit(i));
            out3=repmat(v3(:,i)+ud,1,floor(T3bit(i)*fraction3));
            out4=repmat(v4(:,i)+ud,1,T4bit(i));
            out7=repmat(v3(:,i)+ud,1,T3bit(i)-floor(T3bit(i)*fraction3));
            out012343=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
            out(:,bit*(i-1)+1:bit*i)=out012343;
        case 3%743212
            out0=repmat([1;1;1;1;1]+ud,1,Tzbit(i));
            out1=repmat(v4(:,i)+ud,1,T4bit(i));
            out2=repmat(v3(:,i)+ud,1,T3bit(i));
            out3=repmat(v2(:,i)+ud,1,floor(T2bit(i)*fraction3));
            out4=repmat(v1(:,i)+ud,1,T1bit(i));
            out7=repmat(v2(:,i)+ud,1,T2bit(i)-floor(T2bit(i)*fraction3));
            out743212=[out0 out1 out2 out3 out4 out7 out7 out4 out3 out2 out1 out0];
            out(:,bit*(i-1)+1:bit*i)=out743212;
   end
end
out_2 = circshift(out,[0 1]);
out3 = [out;out_2];
diff1 = sum(abs(out(1,:)-out_2(1,:)));
diff2 = sum(abs(out(2,:)-out_2(2,:)));
diff3 = sum(abs(out(3,:)-out_2(3,:)));
fprintf('out1�������� = %f\n',diff1);
fprintf('out2�������� = %f\n',diff2);
fprintf('out3�������� = %f\n',diff3);

N = time*fsw*bit;
Van=[4/5 -1/5 -1/5 -1/5 -1/5]*out;
Vbn=[-1/5 4/5 -1/5 -1/5 -1/5]*out;
Vcn=[-1/5 -1/5 4/5 -1/5 -1/5]*out;
Vdn=[-1/5 -1/5 -1/5 4/5 -1/5]*out;
Ven=[-1/5 -1/5 -1/5 -1/5 4/5]*out;
fft_Van=fft(Van,N)/N;%6:0.499�B��ݳ�����
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
%i=6:5:�줤��
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(���W^2/���W^2)
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
%ITHD=�W�v���H������
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%������H���t+1=[(N/2-4)-6]/5+1
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
%�Ӻ⪺,����H�W����X�ӤF,���U�ӬO�e��------------------------------------
n2=2:1:N/2;%�h�Y(�Y�Ĥ@�ӭȦ��~���s)
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
disperence=[point3;allarea3;THD_Van;ITHD_Van];
ITHDavg = (ITHD_Van+ITHD_Vbn+ITHD_Vcn)/3;

figure;
subplot(6,1,1),plot(t,s1,t,s2,t,s3,t,s4,t,s5);
subplot(6,1,2),plot(out(1,:));
subplot(6,1,3),plot(out(2,:));
subplot(6,1,4),plot(out(3,:));
subplot(6,1,5),plot(out(4,:));
subplot(6,1,6),plot(out(5,:));

figure;%allHDF
angle=atan2(sindq_p1y,sindq_p1x);
for i=1:1:n
    if angle(i)<0
        angle(i)=angle(i)+2*pi;
    end
end
plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
plot(angle*180/pi,HDF0121,'b'),hold on;
plot(angle*180/pi,HDF7212,'g'),hold on;
% plot(angle*180/pi,minHDFvalue3,'m'),hold on;
%}

% figure,plot(sindqx,sindqy),axis([-0.5 0.5 -0.5 0.5])