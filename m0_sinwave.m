clear
close all
clc
%{
3�Ӭۦ쪺sin�i
t = linspace(0,0.1);				% �b 0 �� 2*pi ���A������ 100 ���I  
a=0.5;
f=50;
bit=300;
fsw=5000;
y=a*sin(2*pi*f*t);
b=a*sin(2*pi*f*t+2*pi/3);
c=a*sin(2*pi*f*t+4*pi/3);
t1=0:1/fsw/bit:0.1;
d=a*sawtooth(2*pi*fsw*t1,0.5);
plot(t, y,t,b,t,c,t1,d)

%}
%{
�P�ɥX�{3�i��
t = linspace(0,0.02);				
a=0.3;
y=a*sin(2*pi*50*t);
b=a*sin(2*pi*100*t);
c=a*sin(2*pi*400*t);
subplot(3,1,1),plot(t,y);
subplot(3,1,2),plot(t,b);
subplot(3,1,3),plot(t,c);
%}

%���X3��sin�i�üе��O
t = 0:1/50000:0.02;				
a=0.3;
y=a*sin(2*pi*50*t);
b=a*sin(2*pi*100*t);
c=a*sin(2*pi*400*t);
plot(t,y,'r:',t,b,'g--',t,c,'b-');
legend('0.3*sin(2*pi*50*t)','0.3*sin(2*pi*100*t)','0.3*sin(2*pi*400*t)'); 
text(0.015,0.3*sin((2*pi*50*0.015)),'\leftarrow0.3sin(100*pi*t)');
text(0.0125,0.3*sin(2*pi*100*0.0125),'\leftarrow0.3sin(200*pi*t)');
text(0.004,0.3*sin((2*pi*400*0.004)),'0.3sin(200*pi*t)\rightarrow'); 








