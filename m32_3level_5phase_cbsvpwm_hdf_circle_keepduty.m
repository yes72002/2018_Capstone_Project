clear%4/29五相
close all
clc
start_time = clock;
time = 0.02;
fsw = 15000;
bit = 400;
amp = 2;%0.2472,0.4,0.6472
fraction3 = 0.5;
f = 50;
n = fsw*time;
t = 0:1/fsw:time-1/fsw;
interval = 0.01;
r = fix(amp/interval);
for j=1:1:r
    amp = interval*(j-1);
    s1(j,:) = amp*cos(2*pi*f*t);
    s2(j,:) = amp*cos(2*pi*f*t-2*pi/5); 
    s3(j,:) = amp*cos(2*pi*f*t-4*pi/5);
    s4(j,:) = amp*cos(2*pi*f*t-6*pi/5); 
    s5(j,:) = amp*cos(2*pi*f*t-8*pi/5); 
end
a = 2*pi/5;
mapping=2/5*[...
    cos(0) cos(a) cos(2*a) cos(3*a) cos(4*a);...
    sin(0) sin(a) sin(2*a) sin(3*a) sin(4*a);...
    cos(0) cos(2*a) cos(4*a) cos(6*a) cos(8*a);...
    sin(0) sin(2*a) sin(4*a) sin(6*a) sin(8*a);...
    1 1 1 1 1];
for j=1:1:r
sindq(5*j-4,:) = mapping(1,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j-3,:) = mapping(2,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j-2,:) = mapping(3,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j-1,:) = mapping(4,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq(5*j,:) = mapping(5,:)*[s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sindq_p1x(j,:) = sindq(5*j-4,:);
sindq_p1y(j,:) = sindq(5*j-3,:);
sindq_p2x(j,:) = sindq(5*j-2,:);
sindq_p2y(j,:) = sindq(5*j-1,:);
end
for j=1:1:r
allsin = [s1(j,:);s2(j,:);s3(j,:);s4(j,:);s5(j,:)];
sinmax(j,:) = max(allsin);
sinmin(j,:) = min(allsin);

w=(-max(allsin)-min(allsin))/2;
allsin=allsin+[w;w;w;w;w];
duty = allsin-floor(allsin);

[dutysort,dutycase] = sort(duty);
duty5th(1,:) = dutysort(1,:);
duty4th(1,:) = dutysort(2,:);
duty3rd(1,:) = dutysort(3,:);
duty2nd(1,:) = dutysort(4,:);
duty1st(1,:) = dutysort(5,:);
Duty0(j,:) = 1-duty1st;
Duty1(j,:) = duty1st-duty2nd;
Duty2(j,:) = duty2nd-duty3rd;
Duty3(j,:) = duty3rd-duty4th;
Duty4(j,:) = duty4th-duty5th;
Duty5(j,:) = duty5th;
updown = floor(allsin);
for i = 1:fsw*time
    v(dutycase(1,i),:) = [0 0 0 0 0 1];
    v(dutycase(2,i),:) = [0 0 0 0 1 1];
    v(dutycase(3,i),:) = [0 0 0 1 1 1];
    v(dutycase(4,i),:) = [0 0 1 1 1 1];
    v(dutycase(5,i),:) = [0 1 1 1 1 1];
    v0(5*j-4:5*j,i) = v(:,1)+updown(:,i);
    v1(5*j-4:5*j,i) = v(:,2)+updown(:,i);
    v2(5*j-4:5*j,i) = v(:,3)+updown(:,i);
    v3(5*j-4:5*j,i) = v(:,4)+updown(:,i);
    v4(5*j-4:5*j,i) = v(:,5)+updown(:,i);
    v5(5*j-4:5*j,i) = v(:,6)+updown(:,i);
    A = [4/5 -1/5 -1/5 -1/5 -1/5;...
    -1/5 4/5 -1/5 -1/5 -1/5;...
    -1/5 -1/5 4/5 -1/5 -1/5;...
    -1/5 -1/5 -1/5 4/5 -1/5;...
    -1/5 -1/5 -1/5 -1/5 4/5];
    vdq0(5*j-4:5*j,i) = mapping*A*v0(5*j-4:5*j,i);
    vdq1(5*j-4:5*j,i) = mapping*A*v1(5*j-4:5*j,i);
    vdq2(5*j-4:5*j,i) = mapping*A*v2(5*j-4:5*j,i);
    vdq3(5*j-4:5*j,i) = mapping*A*v3(5*j-4:5*j,i);
    vdq4(5*j-4:5*j,i) = mapping*A*v4(5*j-4:5*j,i);
    vdq5(5*j-4:5*j,i) = mapping*A*v5(5*j-4:5*j,i);
end
end
for j=1:1:r
    vdq0_p1x(j,:) = vdq0(5*j-4,:);
    vdq0_p1y(j,:) = vdq0(5*j-3,:);
    vdq0_p2x(j,:) = vdq0(5*j-2,:);
    vdq0_p2y(j,:) = vdq0(5*j-1,:);
    vdq1_p1x(j,:) = vdq1(5*j-4,:);
    vdq1_p1y(j,:) = vdq1(5*j-3,:);
    vdq1_p2x(j,:) = vdq1(5*j-2,:);
    vdq1_p2y(j,:) = vdq1(5*j-1,:);
    vdq2_p1x(j,:) = vdq2(5*j-4,:);
    vdq2_p1y(j,:) = vdq2(5*j-3,:);
    vdq2_p2x(j,:) = vdq2(5*j-2,:);
    vdq2_p2y(j,:) = vdq2(5*j-1,:);
    vdq3_p1x(j,:) = vdq3(5*j-4,:);
    vdq3_p1y(j,:) = vdq3(5*j-3,:);
    vdq3_p2x(j,:) = vdq3(5*j-2,:);
    vdq3_p2y(j,:) = vdq3(5*j-1,:);
    vdq4_p1x(j,:) = vdq4(5*j-4,:);
    vdq4_p1y(j,:) = vdq4(5*j-3,:);
    vdq4_p2x(j,:) = vdq4(5*j-2,:);
    vdq4_p2y(j,:) = vdq4(5*j-1,:);
    vdq5_p1x(j,:) = vdq5(5*j-4,:);
    vdq5_p1y(j,:) = vdq5(5*j-3,:);
    vdq5_p2x(j,:) = vdq5(5*j-2,:);
    vdq5_p2y(j,:) = vdq5(5*j-1,:);
end
for j=1:r
    X0_p1(j,:) = (vdq0_p1x(j,:)-sindq_p1x(j,:)).*Duty0(j,:);%v0在dq平面1上x軸
    X1_p1(j,:) = (vdq1_p1x(j,:)-sindq_p1x(j,:)).*Duty1(j,:);%v1在dq平面1上x軸
    X2_p1(j,:) = (vdq2_p1x(j,:)-sindq_p1x(j,:)).*Duty2(j,:);%v2在dq平面1上x軸
    X3_p1(j,:) = (vdq3_p1x(j,:)-sindq_p1x(j,:)).*Duty3(j,:);%v3在dq平面1上x軸
    X4_p1(j,:) = (vdq4_p1x(j,:)-sindq_p1x(j,:)).*Duty4(j,:);%v4在dq平面1上x軸
    X5_p1(j,:) = (vdq5_p1x(j,:)-sindq_p1x(j,:)).*Duty5(j,:);%v5在dq平面1上x軸
    Y0_p1(j,:) = (vdq0_p1y(j,:)-sindq_p1y(j,:)).*Duty0(j,:);%v0在dq平面1上y軸
    Y1_p1(j,:) = (vdq1_p1y(j,:)-sindq_p1y(j,:)).*Duty1(j,:);%v1在dq平面1上y軸
    Y2_p1(j,:) = (vdq2_p1y(j,:)-sindq_p1y(j,:)).*Duty2(j,:);%v2在dq平面1上y軸
    Y3_p1(j,:) = (vdq3_p1y(j,:)-sindq_p1y(j,:)).*Duty3(j,:);%v3在dq平面1上y軸
    Y4_p1(j,:) = (vdq4_p1y(j,:)-sindq_p1y(j,:)).*Duty4(j,:);%v4在dq平面1上y軸
    Y5_p1(j,:) = (vdq5_p1y(j,:)-sindq_p1y(j,:)).*Duty5(j,:);%v5在dq平面1上y軸
    X0_p2(j,:) = (vdq0_p2x(j,:)-sindq_p2x(j,:)).*Duty0(j,:);%v0在dq平面2上x軸
    X1_p2(j,:) = (vdq1_p2x(j,:)-sindq_p2x(j,:)).*Duty1(j,:);%v1在dq平面2上x軸
    X2_p2(j,:) = (vdq2_p2x(j,:)-sindq_p2x(j,:)).*Duty2(j,:);%v2在dq平面2上x軸
    X3_p2(j,:) = (vdq3_p2x(j,:)-sindq_p2x(j,:)).*Duty3(j,:);%v3在dq平面2上x軸
    X4_p2(j,:) = (vdq4_p2x(j,:)-sindq_p2x(j,:)).*Duty4(j,:);%v4在dq平面2上x軸
    X5_p2(j,:) = (vdq5_p2x(j,:)-sindq_p2x(j,:)).*Duty5(j,:);%v5在dq平面2上x軸
    Y0_p2(j,:) = (vdq0_p2y(j,:)-sindq_p2y(j,:)).*Duty0(j,:);%v0在dq平面2上y軸
    Y1_p2(j,:) = (vdq1_p2y(j,:)-sindq_p2y(j,:)).*Duty1(j,:);%v1在dq平面2上y軸
    Y2_p2(j,:) = (vdq2_p2y(j,:)-sindq_p2y(j,:)).*Duty2(j,:);%v2在dq平面2上y軸
    Y3_p2(j,:) = (vdq3_p2y(j,:)-sindq_p2y(j,:)).*Duty3(j,:);%v3在dq平面2上y軸
    Y4_p2(j,:) = (vdq4_p2y(j,:)-sindq_p2y(j,:)).*Duty4(j,:);%v4在dq平面2上y軸
    Y5_p2(j,:) = (vdq5_p2y(j,:)-sindq_p2y(j,:)).*Duty5(j,:);%v5在dq平面2上y軸
end
Xz_p1 = X0_p1 + X5_p1;
Xz_p2 = X0_p2 + X5_p2;
Yz_p1 = Y0_p1 + Y5_p1;
Yz_p2 = Y0_p2 + Y5_p2;
T0 = Duty0;
T1 = Duty1;%無限跟有限精準度
T2 = Duty2;
T3 = Duty3;
T4 = Duty4;
T5 = Duty5;
Tz = Duty0 + Duty5;
for j=1:r
for i=1:fsw*time%Plane1
    P1 = 0;%HDF012347
    P2 = X0_p1(j,i);
    P3 = X0_p1(j,i)+X1_p1(j,i);
    P4 = X0_p1(j,i)+X1_p1(j,i)+X2_p1(j,i);
    P5 = X0_p1(j,i)+X1_p1(j,i)+X2_p1(j,i)+X3_p1(j,i);
    P6 = X0_p1(j,i)+X1_p1(j,i)+X2_p1(j,i)+X3_p1(j,i)+X4_p1(j,i);
    P7 = X0_p1(j,i)+X1_p1(j,i)+X2_p1(j,i)+X3_p1(j,i)+X4_p1(j,i)+X5_p1(j,i);
    HDFX_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*T0(j,i)...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i);
    R1 = 0;
    R2 = Y0_p1(j,i);
    R3 = Y0_p1(j,i)+Y1_p1(j,i);
    R4 = Y0_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i);
    R5 = Y0_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i)+Y3_p1(j,i);
    R6 = Y0_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i)+Y3_p1(j,i)+Y4_p1(j,i);
    R7 = Y0_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i)+Y3_p1(j,i)+Y4_p1(j,i)+Y5_p1(j,i);
    HDFY_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*T0(j,i)...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i);
    P1 = 0;%HDF012343
    P2 = Xz_p1(j,i);
    P3 = Xz_p1(j,i)+X1_p1(j,i);
    P4 = Xz_p1(j,i)+X1_p1(j,i)+X2_p1(j,i);
    P5 = Xz_p1(j,i)+X1_p1(j,i)+X2_p1(j,i)+X3_p1(j,i)*fraction3;
    P6 = Xz_p1(j,i)+X1_p1(j,i)+X2_p1(j,i)+X3_p1(j,i)*fraction3+X4_p1(j,i);
    P7 = Xz_p1(j,i)+X1_p1(j,i)+X2_p1(j,i)+X3_p1(j,i)+X4_p1(j,i);
    HDF0121X_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p1(j,i);
    R3 = Yz_p1(j,i)+Y1_p1(j,i);
    R4 = Yz_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i);
    R5 = Yz_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i)+Y3_p1(j,i)*fraction3;
    R6 = Yz_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i)+Y3_p1(j,i)*fraction3+Y4_p1(j,i);
    R7 = Yz_p1(j,i)+Y1_p1(j,i)+Y2_p1(j,i)+Y3_p1(j,i)+Y4_p1(j,i);
    HDF0121Y_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(j,i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Xz_p1(j,i);
    P3 = Xz_p1(j,i)+X4_p1(j,i);
    P4 = Xz_p1(j,i)+X4_p1(j,i)+X3_p1(j,i);
    P5 = Xz_p1(j,i)+X4_p1(j,i)+X3_p1(j,i)+X2_p1(j,i)*fraction3;
    P6 = Xz_p1(j,i)+X4_p1(j,i)+X3_p1(j,i)+X2_p1(j,i)*fraction3+X1_p1(j,i);
    P7 = Xz_p1(j,i)+X4_p1(j,i)+X3_p1(j,i)+X2_p1(j,i)+X1_p1(j,i);
    HDF7212X_p1(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p1(j,i);
    R3 = Yz_p1(j,i)+Y4_p1(j,i);
    R4 = Yz_p1(j,i)+Y4_p1(j,i)+Y3_p1(j,i);
    R5 = Yz_p1(j,i)+Y4_p1(j,i)+Y3_p1(j,i)+Y2_p1(j,i)*fraction3;
    R6 = Yz_p1(j,i)+Y4_p1(j,i)+Y3_p1(j,i)+Y2_p1(j,i)*fraction3+Y1_p1(j,i);
    R7 = Yz_p1(j,i)+Y4_p1(j,i)+Y3_p1(j,i)+Y2_p1(j,i)+Y1_p1(j,i);
    HDF7212Y_p1(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*(1-fraction3);
end
end
for j=1:r
for i=1:fsw*time%Plane2
    P1 = 0;%HDF012347
    P2 = X0_p2(j,i);
    P3 = X0_p2(j,i)+X1_p2(j,i);
    P4 = X0_p2(j,i)+X1_p2(j,i)+X2_p2(j,i);
    P5 = X0_p2(j,i)+X1_p2(j,i)+X2_p2(j,i)+X3_p2(j,i);
    P6 = X0_p2(j,i)+X1_p2(j,i)+X2_p2(j,i)+X3_p2(j,i)+X4_p2(j,i);
    P7 = X0_p2(j,i)+X1_p2(j,i)+X2_p2(j,i)+X3_p2(j,i)+X4_p2(j,i)+X5_p2(j,i);
    HDFX_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*T0(j,i) ...
              +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
              +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
              +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)...
              +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
              +((P6)^2+(P6)*(P7)+(P7)^2)*T5(j,i);
    R1 = 0;
    R2 = Y0_p2(j,i);
    R3 = Y0_p2(j,i)+Y1_p2(j,i);
    R4 = Y0_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i);
    R5 = Y0_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i)+Y3_p2(j,i);
    R6 = Y0_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i)+Y3_p2(j,i)+Y4_p2(j,i);
    R7 = Y0_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i)+Y3_p2(j,i)+Y4_p2(j,i)+Y5_p2(j,i);
    HDFY_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*T0(j,i) ...
              +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
              +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
              +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)...
              +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
              +((R6)^2+(R6)*(R7)+(R7)^2)*T5(j,i);
    P1 = 0;%HDF012343
    P2 = Xz_p2(j,i);
    P3 = Xz_p2(j,i)+X1_p2(j,i);
    P4 = Xz_p2(j,i)+X1_p2(j,i)+X2_p2(j,i);
    P5 = Xz_p2(j,i)+X1_p2(j,i)+X2_p2(j,i)+X3_p2(j,i)*fraction3;
    P6 = Xz_p2(j,i)+X1_p2(j,i)+X2_p2(j,i)+X3_p2(j,i)*fraction3+X4_p2(j,i);
    P7 = Xz_p2(j,i)+X1_p2(j,i)+X2_p2(j,i)+X3_p2(j,i)+X4_p2(j,i);
    HDF0121X_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T4(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T3(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p2(j,i);
    R3 = Yz_p2(j,i)+Y1_p2(j,i);
    R4 = Yz_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i);
    R5 = Yz_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i)+Y3_p2(j,i)*fraction3;
    R6 = Yz_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i)+Y3_p2(j,i)*fraction3+Y4_p2(j,i);
    R7 = Yz_p2(j,i)+Y1_p2(j,i)+Y2_p2(j,i)+Y3_p2(j,i)+Y4_p2(j,i);
    HDF0121Y_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T4(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T3(j,i)*(1-fraction3);
    P1 = 0;%HDF743212
    P2 = Xz_p2(j,i);
    P3 = Xz_p2(j,i)+X4_p2(j,i);
    P4 = Xz_p2(j,i)+X4_p2(j,i)+X3_p2(j,i);
    P5 = Xz_p2(j,i)+X4_p2(j,i)+X3_p2(j,i)+X2_p2(j,i)*fraction3;
    P6 = Xz_p2(j,i)+X4_p2(j,i)+X3_p2(j,i)+X2_p2(j,i)*fraction3+X1_p2(j,i);
    P7 = Xz_p2(j,i)+X4_p2(j,i)+X3_p2(j,i)+X2_p2(j,i)+X1_p2(j,i);
    HDF7212X_p2(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
                  +((P2)^2+(P2)*(P3)+(P3)^2)*T4(j,i)...
                  +((P3)^2+(P3)*(P4)+(P4)^2)*T3(j,i)...
                  +((P4)^2+(P4)*(P5)+(P5)^2)*T2(j,i)*fraction3...
                  +((P5)^2+(P5)*(P6)+(P6)^2)*T1(j,i)...
                  +((P6)^2+(P6)*(P7)+(P7)^2)*T2(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Yz_p2(j,i);
    R3 = Yz_p2(j,i)+Y4_p2(j,i);
    R4 = Yz_p2(j,i)+Y4_p2(j,i)+Y3_p2(j,i);
    R5 = Yz_p2(j,i)+Y4_p2(j,i)+Y3_p2(j,i)+Y2_p2(j,i)*fraction3;
    R6 = Yz_p2(j,i)+Y4_p2(j,i)+Y3_p2(j,i)+Y2_p2(j,i)*fraction3+Y1_p2(j,i);
    R7 = Yz_p2(j,i)+Y4_p2(j,i)+Y3_p2(j,i)+Y2_p2(j,i)+Y1_p2(j,i);
    HDF7212Y_p2(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
                  +((R2)^2+(R2)*(R3)+(R3)^2)*T4(j,i)...
                  +((R3)^2+(R3)*(R4)+(R4)^2)*T3(j,i)...
                  +((R4)^2+(R4)*(R5)+(R5)^2)*T2(j,i)*fraction3...
                  +((R5)^2+(R5)*(R6)+(R6)^2)*T1(j,i)...
                  +((R6)^2+(R6)*(R7)+(R7)^2)*T2(j,i)*(1-fraction3);
end
end
HDF = HDFX_p1 + HDFY_p1 + HDFX_p2 + HDFY_p2;
HDF0121 = HDF0121X_p1 + HDF0121Y_p1 + HDF0121X_p2 + HDF0121Y_p2;
HDF7212 = HDF7212X_p1 + HDF7212Y_p1 + HDF7212X_p2 + HDF7212Y_p2;
for j=1:1:r
    threeHDF(3*j-2,:)=HDF(j,:);
    threeHDF(3*j-1,:)=HDF0121(j,:);
    threeHDF(3*j,:)=HDF7212(j,:);
    minHDFvalue(j,:)=min(threeHDF(3*j-2:3*j,:));
end
for j=1:1:r
    for i=1:1:fsw*time
        if threeHDF(3*j-2,i)==minHDFvalue(j,i)
            minHDFcase(j,i)=1;
        elseif threeHDF(3*j-1,i)==minHDFvalue(j,i)
            minHDFcase(j,i)=2;
        elseif threeHDF(3*j,i)==minHDFvalue(j,i)
            minHDFcase(j,i)=3;
        end
    end
end
figure;
for j=1:1:r
    plot(sindq_p1x(j,:),sindq_p1y(j,:),'b.');
    hold on;
end
 title('sinwave');
 axis([-inf inf -inf inf],'square')
figure;%HDF圓的分布
sinHDFdistribute1x=zeros(r,fsw*time);
sinHDFdistribute1y=zeros(r,fsw*time);
sinHDFdistribute2x=zeros(r,fsw*time);
sinHDFdistribute2y=zeros(r,fsw*time);
sinHDFdistribute3x=zeros(r,fsw*time);
sinHDFdistribute3y=zeros(r,fsw*time);
for j=1:1:r
for i=1:1:fsw*time
    switch minHDFcase(j,i)
        case 1
        sinHDFdistribute1x(j,i)=sindq_p1x(j,i);
        sinHDFdistribute1y(j,i)=sindq_p1y(j,i);
        case 2
        sinHDFdistribute2x(j,i)=sindq_p1x(j,i);
        sinHDFdistribute2y(j,i)=sindq_p1y(j,i);
        case 3
        sinHDFdistribute3x(j,i)=sindq_p1x(j,i);
        sinHDFdistribute3y(j,i)=sindq_p1y(j,i);
    end
end
end
for j=1:1:r
plot(sinHDFdistribute1x(j,:),sinHDFdistribute1y(j,:),'r.');
hold on;
plot(sinHDFdistribute2x(j,:),sinHDFdistribute2y(j,:),'b.');
hold on;
plot(sinHDFdistribute3x(j,:),sinHDFdistribute3y(j,:),'k.');
hold on;
axis([-inf inf -inf inf],'square')
end
title('hdf circle');

Fdist=sqrt(sum(HDF'));
Fdist0121=sqrt(sum(HDF0121'));
Fdist7212=sqrt(sum(HDF7212'));
figure,plot(Fdist,'r');
hold on,plot(Fdist0121,'b');
hold on,plot(Fdist7212,'g');
title('Fdist');

end_time=clock;
execution_time=end_time-start_time;