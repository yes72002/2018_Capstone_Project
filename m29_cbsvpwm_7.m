clear%1/28�C��cbsvpwm-��m20�@��
close all
clc
start_time=clock;
time = 0.1;
fsw = 36000;
bit = 400;
amp = 0.48;
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

T1=Duty1;%�L���򦳭���ǫ�
T2=Duty2;
T3=Duty3;
T4=Duty4;
T5=Duty5;
T6=Duty6;
T1bit=round(T1*bit/2);
T2bit=round(T2*bit/2);
T3bit=round(T3*bit/2);
T4bit=round(T4*bit/2);
T5bit=round(T5*bit/2);
T6bit=round(T6*bit/2);
Tzbit=bit/2-T1bit-T2bit-T3bit-T4bit-T5bit-T6bit;
for i=1:fsw*time%�D�n���j�{���X
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
fprintf('out1��������=%f\n',diff1);
fprintf('out2��������=%f\n',diff2);
fprintf('out3��������=%f\n',diff3);
fprintf('out4��������=%f\n',diff4);
fprintf('out5��������=%f\n',diff5);
fprintf('out6��������=%f\n',diff6);
fprintf('out7��������=%f\n',diff7);

N = time*fsw*bit;
Van=[6/7 -1/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vbn=[-1/7 6/7 -1/7 -1/7 -1/7 -1/7 -1/7]*out;
Vcn=[-1/7 -1/7 6/7 -1/7 -1/7 -1/7 -1/7]*out;
Vdn=[-1/7 -1/7 -1/7 6/7 -1/7 -1/7 -1/7]*out;
Ven=[-1/7 -1/7 -1/7 -1/7 6/7 -1/7 -1/7]*out;
Vfn=[-1/7 -1/7 -1/7 -1/7 -1/7 6/7 -1/7]*out;
Vgn=[-1/7 -1/7 -1/7 -1/7 -1/7 -1/7 6/7]*out;
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
fft_Vfn=fft(Vfn,N)/N;
mag_Vfn=abs(fft_Vfn)*2;
fft_Vgn=fft(Vgn,N)/N;
mag_Vgn=abs(fft_Vgn)*2;
k=f/(bit*fsw/N)+1:f/(bit*fsw/N):N/2;
%i=(50/10)+1:(50/10):100000
%i=6:5:�줤��
mag_Van_1=mag_Van(1,k);
mag_Vbn_1=mag_Vbn(1,k);
mag_Vcn_1=mag_Vcn(1,k);
mag_Vdn_1=mag_Vdn(1,k);
mag_Ven_1=mag_Ven(1,k);
mag_Vfn_1=mag_Vfn(1,k);
mag_Vgn_1=mag_Vgn(1,k);
% figure;==================================================================
% semilogx(mag_Van_1,'r'),title('mag Van 1');
%THD=100*sqrt(���W^2/���W^2)
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
%ITHD=�W�v���H������
mag_ithd_out1=zeros(1,N/5/2-1);
mag_ithd_out2=zeros(1,N/5/2-1);
mag_ithd_out3=zeros(1,N/5/2-1);
mag_ithd_out4=zeros(1,N/5/2-1);
mag_ithd_out5=zeros(1,N/5/2-1);
mag_ithd_out6=zeros(1,N/5/2-1);
mag_ithd_out7=zeros(1,N/5/2-1);
%fsw*bit*0.02/2-1%������H���t+1=[(N/2-4)-6]/5+1
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
%�Ӻ⪺,����H�W����X�ӤF,���U�ӬO�e��
n2=2:1:N/2;%�h�Y(�Y�Ĥ@�ӭȦ��~���s)
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