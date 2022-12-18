clear
close all
clc
% run('m8_svpwm_3.m');minHDFcase3=zeros(1,fsw*time);
run('m13_svpwm_3_hdf.m');
% run('m22_svpwm_3_newarea3.m');
% run('m22_svpwm_3_newarea5.m');
svpwm=out;
v1partion=v1;
v2partion=v2;
minHDFcasepartion=minHDFcase3;
% run('m25_cbsvpwm_3.m');minHDFcase3=zeros(1,fsw*time);
run('m25_cbsvpwm_3_hdf.m');
% run('m25_cbsvpwm_3_newarea3.m');
% run('m25_cbsvpwm_3_newarea5.m');
cbsvpwm=out;
v1duty=v1;
v2duty=v2;
minHDFcaseduty=minHDFcase3;
sum(sum(abs(v1partion-v1duty)))
sum(sum(abs(v2partion-v2duty)))
sum(sum(abs(minHDFcasepartion-minHDFcaseduty)))
find(abs(minHDFcasepartion-minHDFcaseduty)==1)
%time=0.1,fsw=36000,bit=0.5,amp=0.5
%0127不一樣的點
%v1=8
%v2=4
%minHDFcase=0
%out=0
%舊三區不一樣的點(fraction3=0.5)
%v1=8
%v2=4
%minHDFcase=14(181,301,421,541,661,781,1021,1741,1861,1981,2101,2341,2461,3181)
%out=728
%新三區不一樣的點(fraction3=0.5)
%v1=8
%v2=4
%minHDFcase=14(181,301,421,541,661,781,1021,1741,1861,1981,2101,2341,2461,3181)
%out=728
%新五區不一樣的點
%v1=8
%v2=4
%minHDFcase=14(181,301,421,541,66,781,1021,1741,1861,1981,2101,2341,2461,3181)
%out=728
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
