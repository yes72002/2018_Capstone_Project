% clear%1/28cbsvpwm穝き跋-蛤m27妓
% close all
% clc
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
sindq_p1x=sindqxy(1,:);
sindq_p1y=sindqxy(2,:);
sindq_p2x=sindqxy(3,:);
sindq_p2y=sindqxy(4,:);
sindq_p3x=sindqxy(5,:);
sindq_p3y=sindqxy(6,:);
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
X1_p1 = (vdq1(1,:)-sindq_p1x).*Duty1;%v1dqキ1x禸
X2_p1 = (vdq2(1,:)-sindq_p1x).*Duty2;%v2dqキ1x禸
X3_p1 = (vdq3(1,:)-sindq_p1x).*Duty3;%v3dqキ1x禸
X4_p1 = (vdq4(1,:)-sindq_p1x).*Duty4;%v4dqキ1x禸
X5_p1 = (vdq5(1,:)-sindq_p1x).*Duty5;%v5dqキ1x禸
X6_p1 = (vdq6(1,:)-sindq_p1x).*Duty6;%v6dqキ1x禸
Xz_p1 = -sindq_p1x.*Dutyz;
Y1_p1 = (vdq1(2,:)-sindq_p1y).*Duty1;%v1dqキ1y禸
Y2_p1 = (vdq2(2,:)-sindq_p1y).*Duty2;%v2dqキ1y禸
Y3_p1 = (vdq3(2,:)-sindq_p1y).*Duty3;%v3dqキ1y禸
Y4_p1 = (vdq4(2,:)-sindq_p1y).*Duty4;%v4dqキ1y禸
Y5_p1 = (vdq5(2,:)-sindq_p1y).*Duty5;%v5dqキ1y禸
Y6_p1 = (vdq6(2,:)-sindq_p1y).*Duty6;%v6dqキ1y禸
Yz_p1 = -sindq_p1y.*Dutyz;
X1_p2 = (vdq1(3,:)-sindq_p2x).*Duty1;%v1dqキ2x禸
X2_p2 = (vdq2(3,:)-sindq_p2x).*Duty2;%v2dqキ2x禸
X3_p2 = (vdq3(3,:)-sindq_p2x).*Duty3;%v3dqキ2x禸
X4_p2 = (vdq4(3,:)-sindq_p2x).*Duty4;%v4dqキ2x禸
X5_p2 = (vdq5(3,:)-sindq_p2x).*Duty5;%v5dqキ2x禸
X6_p2 = (vdq6(3,:)-sindq_p2x).*Duty6;%v6dqキ2x禸
Xz_p2 = -sindq_p2x.*Dutyz;
Y1_p2 = (vdq1(4,:)-sindq_p2y).*Duty1;%v6dqキ2y禸
Y2_p2 = (vdq2(4,:)-sindq_p2y).*Duty2;%v6dqキ2y禸
Y3_p2 = (vdq3(4,:)-sindq_p2y).*Duty3;%v6dqキ2y禸
Y4_p2 = (vdq4(4,:)-sindq_p2y).*Duty4;%v6dqキ2y禸
Y5_p2 = (vdq5(4,:)-sindq_p2y).*Duty5;%v6dqキ2y禸
Y6_p2 = (vdq6(4,:)-sindq_p2y).*Duty6;%v6dqキ2y禸
Yz_p2 = -sindq_p2y.*Dutyz;
X1_p3 = (vdq1(5,:)-sindq_p3x).*Duty1;%v1dqキ3x禸
X2_p3 = (vdq2(5,:)-sindq_p3x).*Duty2;%v2dqキ3x禸
X3_p3 = (vdq3(5,:)-sindq_p3x).*Duty3;%v3dqキ3x禸
X4_p3 = (vdq4(5,:)-sindq_p3x).*Duty4;%v4dqキ3x禸
X5_p3 = (vdq5(5,:)-sindq_p3x).*Duty5;%v5dqキ3x禸
X6_p3 = (vdq6(5,:)-sindq_p3x).*Duty6;%v6dqキ3x禸
Xz_p3 = -sindq_p3x.*Dutyz;
Y1_p3 = (vdq1(6,:)-sindq_p3y).*Duty1;%v6dqキ3y禸
Y2_p3 = (vdq2(6,:)-sindq_p3y).*Duty2;%v6dqキ3y禸
Y3_p3 = (vdq3(6,:)-sindq_p3y).*Duty3;%v6dqキ3y禸
Y4_p3 = (vdq4(6,:)-sindq_p3y).*Duty4;%v6dqキ3y禸
Y5_p3 = (vdq5(6,:)-sindq_p3y).*Duty5;%v6dqキ3y禸
Y6_p3 = (vdq6(6,:)-sindq_p3y).*Duty6;%v6dqキ3y禸
Yz_p3 = -sindq_p3y.*Dutyz;
T1=Duty1;%礚蛤Τ弘非
T2=Duty2;
T3=Duty3;
T4=Duty4;
T5=Duty5;
T6=Duty6;
Tz=Dutyz;
for i=1:1:fsw*time%Plane1
    P1=0;%01234567
    P2=Xz_p1(i)/2;
    P3=Xz_p1(i)/2+X1_p1(i);
    P4=Xz_p1(i)/2+X1_p1(i)+X2_p1(i);
    P5=Xz_p1(i)/2+X1_p1(i)+X2_p1(i)+X3_p1(i);
    P6=Xz_p1(i)/2+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i);
    P7=Xz_p1(i)/2+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i);
    P8=Xz_p1(i)/2+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)+X6_p1(i);
    P9=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)+X6_p1(i);
    HDFX_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
           +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
           +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
           +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
           +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
           +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)...
           +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
           +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Yz_p1(i)/2;
    R3=Yz_p1(i)/2+Y1_p1(i);
    R4=Yz_p1(i)/2+Y1_p1(i)+Y2_p1(i);
    R5=Yz_p1(i)/2+Y1_p1(i)+Y2_p1(i)+Y3_p1(i);
    R6=Yz_p1(i)/2+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i);
    R7=Yz_p1(i)/2+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i);
    R8=Yz_p1(i)/2+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)+Y6_p1(i);
    R9=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)+Y6_p1(i);
    HDFY_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
           +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
           +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
           +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
           +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
           +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)...
           +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
           +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Xz_p1(i);
    P3=Xz_p1(i)+X1_p1(i);
    P4=Xz_p1(i)+X1_p1(i)+X2_p1(i);
    P5=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i);
    P6=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i);
    P7=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)*fraction3;
    P8=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)*fraction3+X6_p1(i);
    P9=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)+X6_p1(i);
    HDF0121X_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
               +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
               +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3...
               +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
               +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Yz_p1(i);
    R3=Yz_p1(i)+Y1_p1(i);
    R4=Yz_p1(i)+Y1_p1(i)+Y2_p1(i);
    R5=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i);
    R6=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i);
    R7=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)*fraction3;
    R8=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)*fraction3+Y6_p1(i);
    R9=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)+Y6_p1(i);
    HDF0121Y_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
               +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
               +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3...
               +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
               +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Xz_p1(i);
    P3=Xz_p1(i)+X6_p1(i);
    P4=Xz_p1(i)+X6_p1(i)+X5_p1(i);
    P5=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i);
    P6=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i);
    P7=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)*fraction3;
    P8=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)*fraction3+X1_p1(i);
    P9=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)+X1_p1(i);
    HDF7212X_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
               +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
               +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3...
               +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
               +((P8)^2+(P8)*(P9)+(P9)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Yz_p1(i);
    R3=Yz_p1(i)+Y6_p1(i);
    R4=Yz_p1(i)+Y6_p1(i)+Y5_p1(i);
    R5=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i);
    R6=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i);
    R7=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)*fraction3;
    R8=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)*fraction3+Y1_p1(i);
    R9=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)+Y1_p1(i);
    HDF7212Y_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
               +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
               +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3...
               +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
               +((R8)^2+(R8)*(R9)+(R9)^2)*T2(i)*(1-fraction3);
    P1=0;%new01234565
    P2=Xz_p1(i);
    P3=Xz_p1(i)+X1_p1(i);
    P4=Xz_p1(i)+X1_p1(i)+X2_p1(i);
    P5=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i);
    P6=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i);
    P7=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)*fraction5;
    P8=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)*fraction5+X6_p1(i);
    P9=Xz_p1(i)+X1_p1(i)+X2_p1(i)+X3_p1(i)+X4_p1(i)+X5_p1(i)+X6_p1(i);
    HDFnew0121X_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction5...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction5);
    R1=0;
    R2=Yz_p1(i);
    R3=Yz_p1(i)+Y1_p1(i);
    R4=Yz_p1(i)+Y1_p1(i)+Y2_p1(i);
    R5=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i);
    R6=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i);
    R7=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)*fraction5;
    R8=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)*fraction5+Y6_p1(i);
    R9=Yz_p1(i)+Y1_p1(i)+Y2_p1(i)+Y3_p1(i)+Y4_p1(i)+Y5_p1(i)+Y6_p1(i);
    HDFnew0121Y_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction5...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction5);
    P1=0;%new76543212
    P2=Xz_p1(i);
    P3=Xz_p1(i)+X6_p1(i);
    P4=Xz_p1(i)+X6_p1(i)+X5_p1(i);
    P5=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i);
    P6=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i);
    P7=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)*fraction5;
    P8=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)*fraction5+X1_p1(i);
    P9=Xz_p1(i)+X6_p1(i)+X5_p1(i)+X4_p1(i)+X3_p1(i)+X2_p1(i)+X1_p1(i);
    HDFnew7212X_p1(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction5...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T2(i)*(1-fraction5);
    R1=0;
    R2=Yz_p1(i);
    R3=Yz_p1(i)+Y6_p1(i);
    R4=Yz_p1(i)+Y6_p1(i)+Y5_p1(i);
    R5=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i);
    R6=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i);
    R7=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)*fraction5;
    R8=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)*fraction5+Y1_p1(i);
    R9=Yz_p1(i)+Y6_p1(i)+Y5_p1(i)+Y4_p1(i)+Y3_p1(i)+Y2_p1(i)+Y1_p1(i);
    HDFnew7212Y_p1(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction5...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T2(i)*(1-fraction5);
end
for i=1:1:fsw*time%Plane2
    P1=0;%01234567
    P2=Xz_p2(i)/2;
    P3=Xz_p2(i)/2+X1_p2(i);
    P4=Xz_p2(i)/2+X1_p2(i)+X2_p2(i);
    P5=Xz_p2(i)/2+X1_p2(i)+X2_p2(i)+X3_p2(i);
    P6=Xz_p2(i)/2+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i);
    P7=Xz_p2(i)/2+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i);
    P8=Xz_p2(i)/2+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)+X6_p2(i);
    P9=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)+X6_p2(i);
    HDFX_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
           +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
           +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
           +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
           +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
           +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)...
           +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
           +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Yz_p2(i)/2;
    R3=Yz_p2(i)/2+Y1_p2(i);
    R4=Yz_p2(i)/2+Y1_p2(i)+Y2_p2(i);
    R5=Yz_p2(i)/2+Y1_p2(i)+Y2_p2(i)+Y3_p2(i);
    R6=Yz_p2(i)/2+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i);
    R7=Yz_p2(i)/2+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i);
    R8=Yz_p2(i)/2+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)+Y6_p2(i);
    R9=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)+Y6_p2(i);
    HDFY_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
           +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
           +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
           +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
           +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
           +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)...
           +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
           +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Xz_p2(i);
    P3=Xz_p2(i)+X1_p2(i);
    P4=Xz_p2(i)+X1_p2(i)+X2_p2(i);
    P5=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i);
    P6=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i);
    P7=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)*fraction3;
    P8=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)*fraction3+X6_p2(i);
    P9=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)+X6_p2(i);
    HDF0121X_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
               +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
               +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3...
               +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
               +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Yz_p2(i);
    R3=Yz_p2(i)+Y1_p2(i);
    R4=Yz_p2(i)+Y1_p2(i)+Y2_p2(i);
    R5=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i);
    R6=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i);
    R7=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)*fraction3;
    R8=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)*fraction3+Y6_p2(i);
    R9=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)+Y6_p2(i);
    HDF0121Y_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
               +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
               +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3...
               +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
               +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Xz_p2(i);
    P3=Xz_p2(i)+X6_p2(i);
    P4=Xz_p2(i)+X6_p2(i)+X5_p2(i);
    P5=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i);
    P6=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i);
    P7=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)*fraction3;
    P8=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)*fraction3+X1_p2(i);
    P9=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)+X1_p2(i);
    HDF7212X_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
               +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
               +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3...
               +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
               +((P8)^2+(P8)*(P9)+(P9)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Yz_p2(i);
    R3=Yz_p2(i)+Y6_p2(i);
    R4=Yz_p2(i)+Y6_p2(i)+Y5_p2(i);
    R5=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i);
    R6=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i);
    R7=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)*fraction3;
    R8=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)*fraction3+Y1_p2(i);
    R9=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)+Y1_p2(i);
    HDF7212Y_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
               +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
               +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3...
               +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
               +((R8)^2+(R8)*(R9)+(R9)^2)*T2(i)*(1-fraction3);
    P1=0;%new01234565
    P2=Xz_p2(i);
    P3=Xz_p2(i)+X1_p2(i);
    P4=Xz_p2(i)+X1_p2(i)+X2_p2(i);
    P5=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i);
    P6=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i);
    P7=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)*fraction5;
    P8=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)*fraction5+X6_p2(i);
    P9=Xz_p2(i)+X1_p2(i)+X2_p2(i)+X3_p2(i)+X4_p2(i)+X5_p2(i)+X6_p2(i);
    HDFnew0121X_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction5...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction5);
    R1=0;
    R2=Yz_p2(i);
    R3=Yz_p2(i)+Y1_p2(i);
    R4=Yz_p2(i)+Y1_p2(i)+Y2_p2(i);
    R5=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i);
    R6=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i);
    R7=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)*fraction5;
    R8=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)*fraction5+Y6_p2(i);
    R9=Yz_p2(i)+Y1_p2(i)+Y2_p2(i)+Y3_p2(i)+Y4_p2(i)+Y5_p2(i)+Y6_p2(i);
    HDFnew0121Y_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction5...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction5);
    P1=0;%new76543212
    P2=Xz_p2(i);
    P3=Xz_p2(i)+X6_p2(i);
    P4=Xz_p2(i)+X6_p2(i)+X5_p2(i);
    P5=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i);
    P6=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i);
    P7=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)*fraction5;
    P8=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)*fraction5+X1_p2(i);
    P9=Xz_p2(i)+X6_p2(i)+X5_p2(i)+X4_p2(i)+X3_p2(i)+X2_p2(i)+X1_p2(i);
    HDFnew7212X_p2(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction5...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T2(i)*(1-fraction5);
    R1=0;
    R2=Yz_p2(i);
    R3=Yz_p2(i)+Y6_p2(i);
    R4=Yz_p2(i)+Y6_p2(i)+Y5_p2(i);
    R5=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i);
    R6=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i);
    R7=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)*fraction5;
    R8=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)*fraction5+Y1_p2(i);
    R9=Yz_p2(i)+Y6_p2(i)+Y5_p2(i)+Y4_p2(i)+Y3_p2(i)+Y2_p2(i)+Y1_p2(i);
    HDFnew7212Y_p2(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction5...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T2(i)*(1-fraction5);
end
for i=1:1:fsw*time%Plane3
    P1=0;%01234567
    P2=Xz_p3(i)/2;
    P3=Xz_p3(i)/2+X1_p3(i);
    P4=Xz_p3(i)/2+X1_p3(i)+X2_p3(i);
    P5=Xz_p3(i)/2+X1_p3(i)+X2_p3(i)+X3_p3(i);
    P6=Xz_p3(i)/2+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i);
    P7=Xz_p3(i)/2+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i);
    P8=Xz_p3(i)/2+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)+X6_p3(i);
    P9=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)+X6_p3(i);
    HDFX_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)/2 ...
           +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
           +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
           +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
           +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
           +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)...
           +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
           +((P8)^2+(P8)*(P9)+(P9)^2)*Tz(i)/2;
    R1=0;
    R2=Yz_p3(i)/2;
    R3=Yz_p3(i)/2+Y1_p3(i);
    R4=Yz_p3(i)/2+Y1_p3(i)+Y2_p3(i);
    R5=Yz_p3(i)/2+Y1_p3(i)+Y2_p3(i)+Y3_p3(i);
    R6=Yz_p3(i)/2+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i);
    R7=Yz_p3(i)/2+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i);
    R8=Yz_p3(i)/2+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)+Y6_p3(i);
    R9=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)+Y6_p3(i);
    HDFY_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)/2 ...
           +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
           +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
           +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
           +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
           +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)...
           +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
           +((R8)^2+(R8)*(R9)+(R9)^2)*Tz(i)/2;
    P1=0;%01234565
    P2=Xz_p3(i);
    P3=Xz_p3(i)+X1_p3(i);
    P4=Xz_p3(i)+X1_p3(i)+X2_p3(i);
    P5=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i);
    P6=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i);
    P7=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)*fraction3;
    P8=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)*fraction3+X6_p3(i);
    P9=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)+X6_p3(i);
    HDF0121X_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
               +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
               +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction3...
               +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
               +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction3);
    R1=0;
    R2=Yz_p3(i);
    R3=Yz_p3(i)+Y1_p3(i);
    R4=Yz_p3(i)+Y1_p3(i)+Y2_p3(i);
    R5=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i);
    R6=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i);
    R7=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)*fraction3;
    R8=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)*fraction3+Y6_p3(i);
    R9=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)+Y6_p3(i);
    HDF0121Y_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
               +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
               +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction3...
               +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
               +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction3);
    P1=0;%76543212
    P2=Xz_p3(i);
    P3=Xz_p3(i)+X6_p3(i);
    P4=Xz_p3(i)+X6_p3(i)+X5_p3(i);
    P5=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i);
    P6=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i);
    P7=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i)+X2_p3(i)*fraction3;
    P8=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i)+X2_p3(i)*fraction3+X1_p3(i);
    P9=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i)+X2_p3(i)+X1_p3(i);
    HDF7212X_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
               +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
               +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
               +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
               +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
               +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction3...
               +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
               +((P8)^2+(P8)*(P9)+(P9)^2)*T2(i)*(1-fraction3);
    R1=0;
    R2=Yz_p3(i);
    R3=Yz_p3(i)+Y6_p3(i);
    R4=Yz_p3(i)+Y6_p3(i)+Y5_p3(i);
    R5=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i);
    R6=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i);
    R7=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i)+Y2_p3(i)*fraction3;
    R8=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i)+Y2_p3(i)*fraction3+Y1_p3(i);
    R9=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i)+Y2_p3(i)+Y1_p3(i);
    HDF7212Y_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
               +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
               +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
               +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
               +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
               +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction3...
               +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
               +((R8)^2+(R8)*(R9)+(R9)^2)*T2(i)*(1-fraction3);
    P1=0;%new01234565
    P2=Xz_p3(i);
    P3=Xz_p3(i)+X1_p3(i);
    P4=Xz_p3(i)+X1_p3(i)+X2_p3(i);
    P5=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i);
    P6=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i);
    P7=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)*fraction5;
    P8=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)*fraction5+X6_p3(i);
    P9=Xz_p3(i)+X1_p3(i)+X2_p3(i)+X3_p3(i)+X4_p3(i)+X5_p3(i)+X6_p3(i);
    HDFnew0121X_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T5(i)*fraction5...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T6(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T5(i)*(1-fraction5);
    R1=0;
    R2=Yz_p3(i);
    R3=Yz_p3(i)+Y1_p3(i);
    R4=Yz_p3(i)+Y1_p3(i)+Y2_p3(i);
    R5=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i);
    R6=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i);
    R7=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)*fraction5;
    R8=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)*fraction5+Y6_p3(i);
    R9=Yz_p3(i)+Y1_p3(i)+Y2_p3(i)+Y3_p3(i)+Y4_p3(i)+Y5_p3(i)+Y6_p3(i);
    HDFnew0121Y_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T5(i)*fraction5...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T6(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T5(i)*(1-fraction5);
    P1=0;%new76543212
    P2=Xz_p3(i);
    P3=Xz_p3(i)+X6_p3(i);
    P4=Xz_p3(i)+X6_p3(i)+X5_p3(i);
    P5=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i);
    P6=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i);
    P7=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i)+X2_p3(i)*fraction5;
    P8=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i)+X2_p3(i)*fraction5+X1_p3(i);
    P9=Xz_p3(i)+X6_p3(i)+X5_p3(i)+X4_p3(i)+X3_p3(i)+X2_p3(i)+X1_p3(i);
    HDFnew7212X_p3(i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T6(i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T5(i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T4(i)...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T3(i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(i)*fraction5...
                  +((P7)^2+(P7)*(P8)+(P8)^2)*T1(i)...
                  +((P8)^2+(P8)*(P9)+(P9)^2)*T2(i)*(1-fraction5);
    R1=0;
    R2=Yz_p3(i);
    R3=Yz_p3(i)+Y6_p3(i);
    R4=Yz_p3(i)+Y6_p3(i)+Y5_p3(i);
    R5=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i);
    R6=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i);
    R7=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i)+Y2_p3(i)*fraction5;
    R8=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i)+Y2_p3(i)*fraction5+Y1_p3(i);
    R9=Yz_p3(i)+Y6_p3(i)+Y5_p3(i)+Y4_p3(i)+Y3_p3(i)+Y2_p3(i)+Y1_p3(i);
    HDFnew7212Y_p3(i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T6(i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T5(i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T4(i)...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T3(i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(i)*fraction5...
                  +((R7)^2+(R7)*(R8)+(R8)^2)*T1(i)...
                  +((R8)^2+(R8)*(R9)+(R9)^2)*T2(i)*(1-fraction5);
end
HDF = HDFX_p1 + HDFY_p1 + HDFX_p2 + HDFY_p2 + HDFX_p3 + HDFY_p3;
HDF0121 = HDF0121X_p1 + HDF0121Y_p1 + HDF0121X_p2 + HDF0121Y_p2 + HDF0121X_p3 + HDF0121Y_p3;
HDF7212 = HDF7212X_p1 + HDF7212Y_p1 + HDF7212X_p2 + HDF7212Y_p2 + HDF7212X_p3 + HDF7212Y_p3;
HDFnew0121 = HDFnew0121X_p1 + HDFnew0121Y_p1 + HDFnew0121X_p2 + HDFnew0121Y_p2 + HDFnew0121X_p3 + HDFnew0121Y_p3;
HDFnew7212 = HDFnew7212X_p1 + HDFnew7212Y_p1 + HDFnew7212X_p2 + HDFnew7212Y_p2 + HDFnew7212X_p3 + HDFnew7212Y_p3;

threeHDF=[HDF;HDF0121;HDF7212];
fiveHDF=[HDF;HDF0121;HDF7212;HDFnew0121;HDFnew7212];
[minHDFvalue3,minHDFcase3]=min(threeHDF);
[minHDFvalue5,minHDFcase5]=min(fiveHDF);
point3=sum(minHDFcase3==2,2)+sum(minHDFcase3==3,2);
point5=sum(minHDFcase5==4,2)+sum(minHDFcase5==5,2);
allarea3=sum(abs(HDF-minHDFvalue3));
allarea5=sum(abs(minHDFvalue3-minHDFvalue5));

T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
T5bit=round(T5*bit/2);
T6bit=round(T6*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit-T5bit-T6bit;
for i=1:fsw*time%璶祘Α絏
   switch minHDFcase5(i)
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
       case 4%new01234565
           out0=repmat([0;0;0;0;0;0;0],1,Tzbit(i));
           out1=repmat(v1(:,i),1,T1bit(i));
           out2=repmat(v2(:,i),1,T2bit(i));
           out3=repmat(v3(:,i),1,T3bit(i));
           out4=repmat(v4(:,i),1,T4bit(i));
           out5=repmat(v5(:,i),1,floor(T5bit(i)*fraction5));
           out6=repmat(v6(:,i),1,T6bit(i));
           out7=repmat(v5(:,i),1,T5bit(i)-floor(T5bit(i)*fraction5));
           outnew01234565=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=outnew01234565;
       case 5%new76543212
           out0=repmat([1;1;1;1;1;1;1],1,Tzbit(i));
           out1=repmat(v6(:,i),1,T6bit(i));
           out2=repmat(v5(:,i),1,T5bit(i));
           out3=repmat(v4(:,i),1,T4bit(i));
           out4=repmat(v3(:,i),1,T3bit(i));
           out5=repmat(v2(:,i),1,floor(T2bit(i)*fraction5));
           out6=repmat(v1(:,i),1,T1bit(i));
           out7=repmat(v2(:,i),1,T2bit(i)-floor(T2bit(i)*fraction5));
           outnew76543212=[out0 out1 out2 out3 out4 out5 out6 out7 out7 out6 out5 out4 out3 out2 out1 out0];
           out(:,bit*(i-1)+1:bit*i)=outnew76543212;
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
fprintf('out1ち传Ω计=%f\n',diff1);
fprintf('out2ち传Ω计=%f\n',diff2);
fprintf('out3ち传Ω计=%f\n',diff3);
fprintf('out4ち传Ω计=%f\n',diff4);
fprintf('out5ち传Ω计=%f\n',diff5);
fprintf('out6ち传Ω计=%f\n',diff6);
fprintf('out7ち传Ω计=%f\n',diff7);

N = time*fsw*bit;
Van=[6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vbn=[-1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vcn=[-1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7]*out;
Vdn=[-1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7]*out;
Ven=[-1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7]*out;
Vfn=[-1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7]*out;
Vgn=[-1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7]*out;
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
fft_Vfn=fft(Vfn,N)/N;
mag_Vfn=abs(fft_Vfn)*2;
fft_Vgn=fft(Vgn,N)/N;
mag_Vgn=abs(fft_Vgn)*2;
k=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:い丁
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
mag_Vfn_1=mag_Vfn(1,k);
mag_Vgn_1=mag_Vgn(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(繵^2/膀繵^2)
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
%ITHD=繵瞯埃セō
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
mag_ithd_out6=zeros(1,N/5/2-1);
mag_ithd_out7=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%ソ搭埃そ畉+1=[(N/2-4)-6]/5+1
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
%赣衡,︽常衡ㄓ,钡ㄓ琌礶瓜
n2=2:1:N/2;%繷(繷材Τ粇箂)
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
disperence=[point5;allarea5;THD_Van;ITHD_Van;ITHD_Vbn;ITHDavg];
% figure;%allHDF
% angle=atan2(sindq_p1y,sindq_p1x);
% for i=1:1:fsw*time
%     if angle(i)<0
%         angle(i)=angle(i)+2*pi;
%     end
% end
% plot(angle*180/pi,HDF,'r'),title('allHDF'),hold on;
% plot(angle*180/pi,HDF0121,'b'),hold on;
% plot(angle*180/pi,HDF7212,'g'),hold on;
% plot(angle*180/pi,minHDFvalue3,'m'),hold on;
% plot(angle*180/pi,HDFnew0121,'c'),hold on;
% plot(angle*180/pi,HDFnew7212,'y'),hold on;
% plot(angle*180/pi,minHDFvalue5,'k'),hold on;