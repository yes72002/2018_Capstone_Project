clear
close all
clc
%%
% integral---------------int(f,x,�U��,�W��)
% ��x�L��---------------diff(f)
% ��a�L��---------------diff(f, a)
% ��x�G���L��-----------diff(f, 2)
% �]�i��x�}���C�@�Ӥ����i��L��
syms a x;
f=cos(x)^2;
inta=int(f,x,0,2*pi);%pi
g=sin(a*x^2);
dg=diff(g,a);%x^2*cos(a*x^2)
%%
%ln---------------------log
%log--------------------log10
e=log(15.43e11);
f=log10(10^11);
g=log10(100^11);
%%
%�۱�transfer function
Ts=20;
Tp=4;
final=1;
OS=5;
s=-log(OS/100)/sqrt(pi^2+(log(OS/100))^2);
wn=pi/Tp/sqrt(1-s^2);
a=2*s*wn;
b=wn^2;
c=final*wn^2;
s=0.252;
wn=1.789;
Tr=(1.76*s^3-0.417*s^2+1.039*s+1)/wn;
%%
%�h�ۦ��D��
syms s T;
A=[2*s*s/25+2*s/25,-2*s/25,0;-2*s/25,2*s/25+5/25,-5/25;0,-5/25,2*s/25+5/25]*25/9 ;
B=[T;0;0] ;
control=A\B;
%%
%����Ƶ��
%real
%imag
%%
%degree to radius
%deg2rad
%rad2deg
%%
%plot polynomial
%method 1
% y = @(t) t.^-1-t.^-2;
% x = 0:0.1:100;
% figure,plot(x,y(x)),axis([0 10 -1.5 1.5]);
%method 2
ezplot('x^-1-x^-2',[0.1 10]);%�S���q�{-2pi~2pi
axis([0 10 -1.5 1.5])
% method 3%�w�w��
figure,fplot(@(x) x^-1-x^-2,[0.2 10]);
axis([0 10 -1.5 1.5]);
%%
%%�x�}���Ƥ���,�t����0
a=[1 2 -3 4] ;
a.*(a>0)  ;
%%
%�h����
%(4x^4+21x^3+28x^2+8x+3)����F��
poly1=[4 21 28 21 8 3];
%(s-2)(s-7)(s+5)����F��
poly2=poly([2 7 -5]);
%�h�����ۭ�
poly3=conv(poly1,poly2);
%�h�����۰�(polynomial mutiplication)
[poly_q,poly_r]=deconv(poly1,poly2);
%�N�J�h�����᪺�ƭ�
m1=polyval([1 2 3],[1 -1]);
%�h�����L��
poly4=polyder([3 6 9],[2 1]);%(3s^2+6s+7)*(2s+1)���L��
[poly5,poly6]=polyder([-1 -3 -2],[1 -8 15]);%(-s^2-3s-2)/(s^2-8s+15)
%�h�����n��
poly7=polyint([6 6 4],5);%K�O�`��
%�x�}�h�����i�}=>det(sI-A)���Y��
A=[2 1 1;1 7 1;-3 4 -5];
poly(A)
%%
K=0.9;
% numg=[1];
numg=poly([-3 -4]);
% deng1=[1 0  0 0];
deng1=poly([-1 -2]);
% deng1=poly([-10^3 -10^5 -10^6]);
deng2=[1];
G=tf(K*numg,conv(deng1,deng2));
% rlocus(G)
T=G/(1+G);
% ���O-------------------------------------------------------------------------------
% figure,step(T)
% figure,impules(T)
% step=stepinfo(T)
% kp=dcgain(G)
% error_step=1/(kp+1)
% s = tf('s');
% figure,step(T/s);%ramp response
% zpk=[zero,pole,gain];
% rlocus(G)(root locus �S����T,�Х��I��)
% figure,bode(G)
% figure,bode(T)%��resonant frequency��
% [Gm,Pm,Wcg,Wcp]=margin(G)
% nyquist(G)
% axis([-10 10 -10 10])
% A=[2 1 1;1 7 1;-3 4 -5];
% poly(A)%�S�x�h����
% s = tf('s');
% sys = (244.2*s + 244.2) / (0.015*s^4 + 1.525*s^3 + 2.51*s^2 + 245.2*s + 1221);
%Routh-Hurwitz tale
% syms a b c EPS;
% ra=routh([1 a b c],EPS);
% P = tf(num,den,'InputDelay',3.4)%transfer function*exp(-3.4s)
