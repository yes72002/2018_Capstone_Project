clear%三角波比較法
close all
clc
time = 0.1;
fsw = 5000;
bit = 300;
amp = 0.5;
f = 60;
t = 0:1/fsw/bit:time-1/fsw/bit; 
s1 = amp*sin(2*pi*f*t); 
s2 = amp*sin(2*pi*f*t+2*pi/3); 
s3 = amp*sin(2*pi*f*t+4*pi/3); 
tri = amp*sawtooth(2*pi*fsw*t,0.5);
subplot(4,1,1),plot(t,s1,t,s2,t,s3,t,tri);
out1=zeros(1,fsw*bit*time);
for x=1:1:fsw*bit*time
    if s1(x)>=tri(x)
        out1(x)=1;
    else
        out1(x)=0;
    end
end
out2=zeros(1,fsw*bit*time);
for x=1:1:fsw*bit*time
    if s2(x)>=tri(x)
        out2(x)=1;
    else
        out2(x)=0;
    end
end
out3=zeros(1,fsw*bit*time);
for x=1:1:fsw*bit*time
    if s3(x)>=tri(x)
        out3(x)=1;
    else
        out3(x)=0;
    end
end
subplot(4,1,2),plot(0:1:fsw*bit*time-1,out1);
subplot(4,1,3),plot(0:1:fsw*bit*time-1,out2);
subplot(4,1,4),plot(0:1:fsw*bit*time-1,out3);

N=time*fsw*bit;
freqStep = 10;
freq = freqStep*(0:N/2-1);

fft_out1 = fft(out1,N);
mag_out1 = abs(fft_out1)/(length(fft_out1)/2); 
a1=mag_out1(1:1:N/2);

fft_out2 = fft(out2,N);
mag_out2 = abs(fft_out2)/(length(fft_out2)/2); 
a2=mag_out2(1:1:N/2);

fft_out3 = fft(out3,N);
mag_out3 = abs(fft_out3)/(length(fft_out3)/2); 
a3=mag_out3(1:1:N/2);

for i=2:fsw*bit*time/2
  freq(i-1)=freq(i);
  a1(i-1)=a1(i);
  a2(i-1)=a2(i);
  a3(i-1)=a3(i);
end
figure;
subplot(3,1,1),plot(freq,a1)
subplot(3,1,2),plot(freq,a2);
subplot(3,1,3),plot(freq,a3)

harm1 = a1( 12: 6: end/2 ); 
THD1 = 100*sqrt(  sum( harm1.^2) /a1(6)^2 );

harm2 = a2( 12: 6: end/2 ); 
THD2 = 100*sqrt(  sum( harm2.^2) /a2(6)^2 );

harm3 = a3( 12: 6: end/2 ); 
THD3 = 100*sqrt(  sum( harm3.^2) /a3(6)^2 );

A = [ 2/3 -1/3 -1/3 ; -1/3 2/3 -1/3 ; -1/3 -1/3 2/3];  
B = [ out1 ; out2 ; out3 ];  
C = 48*A*B  ; 
figure;
plot(C(1,1:1:fsw*bit*time));
hold on;
plot(C(2,1:1:fsw*bit*time));
hold on;
plot(C(3,1:1:fsw*bit*time));
hold on;
figure
semilogx(C(1,1:1:fsw*bit*time))