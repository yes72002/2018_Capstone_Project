clear
close all
clc
% run('m16_svpwm_5.m');minHDFcase3=zeros(1,fsw*time);
% run('m17_svpwm_5_hdf.m');
run('m23_svpwm_5_newarea3.m');
% run('m23_svpwm_5_newarea5.m');
svpwm=out;
v1partion=v1;
v2partion=v2;
v3partion=v3;
v4partion=v4;
T4partion=T4;
minHDFcasepartion=minHDFcase3;
% run('m27_cbsvpwm_5.m');minHDFcase3=zeros(1,fsw*time);
% % run('m27_cbsvpwm_5_hdf.m');
run('m27_cbsvpwm_5_newarea3.m');
% run('m27_cbsvpwm_5_newarea5.m');
cbsvpwm=out;
v1duty=v1;
v2duty=v2;
v3duty=v3;
v4duty=v4;
minHDFcaseduty=minHDFcase3;
sum(sum(abs(v1partion-v1duty)))
sum(sum(abs(v2partion-v2duty)))
sum(sum(abs(v3partion-v3duty)))
sum(sum(abs(v4partion-v4duty)))
sum(sum(abs(minHDFcasepartion-minHDFcaseduty)))
% find(abs(minHDFcasepartion-minHDFcaseduty)==1)
%time=0.1,fsw=36000,bit=0.5,amp=0.5
%0127不一樣的點
%v1=4,v2=6,v3=2,v4=6
%minHDFcase=0
%out=0
%舊三區不一樣的點
%v1=4,v2=6,v3=2,v4=6
%minHDFcase=0
%out=0
%新三區不一樣的點
%v1=4,v2=6,v3=2,v4=6
%minHDFcase=0
%out=0
%新五區不一樣的點(amp=0.5,fraction3=0.1,fraction5=0.5)
%v1=4,v2=6,v3=2,v4=6
%minHDFcase=0
%out=0
diff=0;
diffout=0;
for i=1:fsw*time*bit
    svout=svpwm(:,i);
    cbout=cbsvpwm(:,i);
    if svout~=cbout
        diff=i;
        diffout=diffout+1;
    end
end
sum(sum(abs(svpwm-cbsvpwm)))
disp(diff)
disp(diffout)
