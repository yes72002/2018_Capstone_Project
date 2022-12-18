clear
close all
clc
% run('m20_svpwm_7.m');minHDFcase3=zeros(1,fsw*time);
% run('m21_svpwm_7_hdf.m');
% run('m24_svpwm_7_newarea3.m');
% run('m24_svpwm_7_newarea5.m');
run('m24_svpwm_7_newarea3_3plane.m');
% run('m24_svpwm_7_newarea5_3plane.m');
svpwm=out;
v1partion=v1;
v2partion=v2;
v3partion=v3;
v4partion=v4;
v5partion=v5;
v6partion=v6;
minHDFcasepartion=minHDFcase3;
% run('m29_cbsvpwm_7.m');minHDFcase3=zeros(1,fsw*time);
% run('m29_cbsvpwm_7_hdf.m');
% run('m29_cbsvpwm_7_newarea3.m');
% run('m29_cbsvpwm_7_newarea5.m');
run('m29_cbsvpwm_7_newarea3_3plane.m');
% run('m29_cbsvpwm_7_newarea5_3plane.m');
cbsvpwm=out;
v1duty=v1;
v2duty=v2;
v3duty=v3;
v4duty=v4;
v5duty=v5;
v6duty=v6;
minHDFcaseduty=minHDFcase3;
disp(sum(sum(abs(v1partion-v1duty))))
disp(sum(sum(abs(v2partion-v2duty))))
disp(sum(sum(abs(v3partion-v3duty))))
disp(sum(sum(abs(v4partion-v4duty))))
disp(sum(sum(abs(v5partion-v5duty))))
disp(sum(sum(abs(v6partion-v6duty))))
sum(sum(abs(minHDFcasepartion-minHDFcaseduty)))
find(abs(minHDFcasepartion-minHDFcaseduty)==1)
%time=0.1,fsw=36000,bit=0.5,amp=0.5
%0127不一樣的點
%v1=2,v2=0,v3=2,v4=0,v5=2,v6=0
%minHDFcase=0
%out=0
%舊三區不一樣的點(amp=0.5)
%v1=2,v2=0,v3=2,v4=0,v5=2,v6=0
%minHDFcase=8(181,541,901,1261,1621,1981,2341,3421)
%out=96
%新三區不一樣的點(amp=0.51)
%v1=2,v2=0,v3=2,v4=0,v5=2,v6=0
%minHDFcase=8(181,901,1261,1621,1981,2341,3061,3421)
%out=32
%新五區不一樣的點
%v1=2,v2=0,v3=2,v4=0,v5=2,v6=0
%minHDFcase=2(901,3421)
%out=8
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
sum(sum(abs(svpwm-cbsvpwm)))%不一樣的點在共幾個點
disp(diff)%不一樣的行在第幾行
disp(diffout)%不一樣的行共幾行
