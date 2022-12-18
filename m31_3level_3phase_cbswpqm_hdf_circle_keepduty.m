clear%3/29三相
close all
clc
start_time = clock;
time  =  0.1;
fsw  =  15000;
bit  =  400;
amp  =  1.44;%0.2472,0.4,0.6472
fraction3 = 0.5;
f  =  50;
n  =  fsw*time;
t  =  0:1/fsw:time-1/fsw;
interval = 0.01;
r = fix(amp/interval);
s1 = zeros(r,fsw*time);
s2 = zeros(r,fsw*time);
s3 = zeros(r,fsw*time);
sindq = zeros(r,fsw*time);
sindqx = zeros(r,fsw*time);
sindqy = zeros(r,fsw*time);
for j=1:1:r
   amp = interval*(j-1);
   s1(j,:) = amp*cos(2*pi*f*t);
   s2(j,:) = amp*cos(2*pi*f*t-2*pi/3);
   s3(j,:) = amp*cos(2*pi*f*t-4*pi/3); 
end
a = 2*pi/3;
mapping = 2/3*[cos(0) cos(a) cos(2*a);sin(0) sin(a) sin(2*a);1 1 1];
for j=1:1:r
sindq(3*j-2,:) = mapping(1,:)*[s1(j,:);s2(j,:);s3(j,:)];
sindq(3*j-1,:) = mapping(2,:)*[s1(j,:);s2(j,:);s3(j,:)];
sindq(3*j,:) = mapping(3,:)*[s1(j,:);s2(j,:);s3(j,:)];
sindqx(j,:) = sindq(3*j-2,:);
sindqy(j,:) = sindq(3*j-1,:);
end
for j=1:1:r
allsin = [s1(j,:);s2(j,:);s3(j,:)];
sinmax(j,:) = max(allsin);
sinmin(j,:) = min(allsin);

w=(-max(allsin)-min(allsin))/2;
allsin=allsin+[w;w;w];
duty = allsin-floor(allsin);

[dutysort,dutycase] = sort(duty);
duty3rd(1,:) = dutysort(1,:);
duty2nd(1,:) = dutysort(2,:);
duty1st(1,:) = dutysort(3,:);
Duty0(j,:) = 1-duty1st;
Duty1(j,:) = duty1st-duty2nd;
Duty2(j,:) = duty2nd-duty3rd;
Duty3(j,:) = duty3rd;
% Dutyz(j,:) = 1-Duty1(j,:)-Duty2(j,:);
updown = floor(allsin);
for i = 1:fsw*time
    v(dutycase(1,i),:) = [0 0 0 1];
    v(dutycase(2,i),:) = [0 0 1 1];
    v(dutycase(3,i),:) = [0 1 1 1];
    v0(3*j-2:3*j,i) = v(:,1)+updown(:,i);
    v1(3*j-2:3*j,i) = v(:,2)+updown(:,i);
    v2(3*j-2:3*j,i) = v(:,3)+updown(:,i);
    v3(3*j-2:3*j,i) = v(:,4)+updown(:,i);
    A = [2/3 -1/3 -1/3;-1/3 2/3 -1/3;-1/3 -1/3 2/3];
    vdq0(3*j-2:3*j,i) = mapping*A*v0(3*j-2:3*j,i);
    vdq1(3*j-2:3*j,i) = mapping*A*v1(3*j-2:3*j,i);
    vdq2(3*j-2:3*j,i) = mapping*A*v2(3*j-2:3*j,i);
    vdq3(3*j-2:3*j,i) = mapping*A*v3(3*j-2:3*j,i);
end
end
for j=1:1:r
    vdq0x(j,:) = vdq0(3*j-2,:);
    vdq0y(j,:) = vdq0(3*j-1,:);
    vdq1x(j,:) = vdq1(3*j-2,:);
    vdq1y(j,:) = vdq1(3*j-1,:);
    vdq2x(j,:) = vdq2(3*j-2,:);
    vdq2y(j,:) = vdq2(3*j-1,:);
    vdq3x(j,:) = vdq3(3*j-2,:);
    vdq3y(j,:) = vdq3(3*j-1,:);
end
vrefx = sindqx;
vrefy = sindqy;

for j=1:r
    X0(j,:) = (vdq0x(j,:)-vrefx(j,:)).*Duty0(j,:);  
    X1(j,:) = (vdq1x(j,:)-vrefx(j,:)).*Duty1(j,:);
    X2(j,:) = (vdq2x(j,:)-vrefx(j,:)).*Duty2(j,:);
    X3(j,:) = (vdq3x(j,:)-vrefx(j,:)).*Duty3(j,:);
    Y0(j,:) = (vdq0y(j,:)-vrefy(j,:)).*Duty0(j,:);
    Y1(j,:) = (vdq1y(j,:)-vrefy(j,:)).*Duty1(j,:);
    Y2(j,:) = (vdq2y(j,:)-vrefy(j,:)).*Duty2(j,:);
    Y3(j,:) = (vdq3y(j,:)-vrefy(j,:)).*Duty3(j,:);
end
Xz = X0 + X3;
Yz = Y0 + Y3;
T0 = Duty0;
T1 = Duty1;%無限跟有限精準度
T2 = Duty2;
T3 = Duty3;
Tz = Duty0 + Duty3;

for j=1:r
for i=1:fsw*time
    P1 = 0;%HDF0127
    P2 = X0(j,i);
    P3 = X0(j,i)+X1(j,i);
    P4 = X0(j,i)+X1(j,i)+X2(j,i);
    P5 = X0(j,i)+X1(j,i)+X2(j,i)+X3(j,i);
    HDFX(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*T0(j,i) ...
          +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)...
          +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
          +((P4)^2+(P4)*(P5)+(P5)^2)*T3(j,i);
    R1 = 0;
    R2 = Y0(j,i);
    R3 = Y0(j,i)+Y1(j,i);
    R4 = Y0(j,i)+Y1(j,i)+Y2(j,i);
    R5 = Y0(j,i)+Y1(j,i)+Y2(j,i)+Y3(j,i);
    HDFY(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*T0(j,i) ...
          +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)...
          +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
          +((R4)^2+(R4)*(R5)+(R5)^2)*T3(j,i);
    P1 = 0;%HDF0121
    P2 = Xz(j,i);
    P3 = Xz(j,i)+X1(j,i)*fraction3;
    P4 = Xz(j,i)+X1(j,i)*fraction3+X2(j,i);
    P5 = Xz(j,i)+X1(j,i)+X2(j,i);
    HDF0121X(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
             +((P2)^2+(P2)*(P3)+(P3)^2)*T1(j,i)*fraction3...
             +((P3)^2+(P3)*(P4)+(P4)^2)*T2(j,i)...
             +((P4)^2+(P4)*(P5)+(P5)^2)*T1(j,i)*(1-fraction3);   
    R1 = 0;
    R2 = Yz(j,i);
    R3 = Yz(j,i)+Y1(j,i)*fraction3;
    R4 = Yz(j,i)+Y1(j,i)*fraction3+Y2(j,i);
    R5 = Yz(j,i)+Y1(j,i)+Y2(j,i);
    HDF0121Y(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
             +((R2)^2+(R2)*(R3)+(R3)^2)*T1(j,i)*fraction3...
             +((R3)^2+(R3)*(R4)+(R4)^2)*T2(j,i)...
             +((R4)^2+(R4)*(R5)+(R5)^2)*T1(j,i)*(1-fraction3);
    P1 = 0;%HDF7212
    P2 = Xz(j,i);
    P3 = Xz(j,i)+X2(j,i)*fraction3;
    P4 = Xz(j,i)+X2(j,i)*fraction3+X1(j,i);
    P5 = Xz(j,i)+X2(j,i)+X1(j,i);
    HDF7212X(j,i)=((P1)^2+(P1)*(P2)+(P2)^2)*Tz(j,i)...
             +((P2)^2+(P2)*(P3)+(P3)^2)*T2(j,i)*fraction3...
             +((P3)^2+(P3)*(P4)+(P4)^2)*T1(j,i)...
             +((P4)^2+(P4)*(P5)+(P5)^2)*T2(j,i)*(1-fraction3);
    R1 = 0;
    R2 = Yz(j,i);
    R3 = Yz(j,i)+Y2(j,i)*fraction3;
    R4 = Yz(j,i)+Y2(j,i)*fraction3+Y1(j,i);
    R5 = Yz(j,i)+Y2(j,i)+Y1(j,i);
    HDF7212Y(j,i)=((R1)^2+(R1)*(R2)+(R2)^2)*Tz(j,i)...
             +((R2)^2+(R2)*(R3)+(R3)^2)*T2(j,i)*fraction3...
             +((R3)^2+(R3)*(R4)+(R4)^2)*T1(j,i)...
             +((R4)^2+(R4)*(R5)+(R5)^2)*T2(j,i)*(1-fraction3);
end
HDF = HDFX + HDFY;
HDF0121 = HDF0121X + HDF0121Y;
HDF7212 = HDF7212X + HDF7212Y;
end
threeHDF=zeros(3*r,fsw*time);
minHDFcase=zeros(r,fsw*time);
minHDFvalue=zeros(r,fsw*time);
for j=1:1:r
    threeHDF(3*j-2,:)=HDF(j,:);
    threeHDF(3*j-1,:)=HDF0121(j,:);
    threeHDF(3*j,:)=HDF7212(j,:);
    minHDFvalue(j,:)=min(threeHDF(3*j-2:3*j,:));
end
figure;
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
for j=1:1:r
    plot(sindqx(j,:),sindqy(j,:),'b.');
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
        sinHDFdistribute1x(j,i)=sindqx(j,i);
        sinHDFdistribute1y(j,i)=sindqy(j,i);
        case 2
        sinHDFdistribute2x(j,i)=sindqx(j,i);
        sinHDFdistribute2y(j,i)=sindqy(j,i);
        case 3
        sinHDFdistribute3x(j,i)=sindqx(j,i);
        sinHDFdistribute3y(j,i)=sindqy(j,i);
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
end
title('hdf circle');
axis([-inf inf -inf inf],'square')

Fdist=sqrt(sum(HDF'));
Fdist0121=sqrt(sum(HDF0121'));
Fdist7212=sqrt(sum(HDF7212'));
figure,plot(Fdist,'r');
hold on,plot(Fdist0121,'b');
hold on,plot(Fdist7212,'g');
title('Fdist');