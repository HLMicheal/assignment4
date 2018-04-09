%ELEC4700 Assignment 4
% Part 2 
% By Huanyu Liu
% 100986552
clear
clc
% parameters
R1=1;
c=0.25;
R2=2;
L=0.2;
R3=10;
alpha=100;
R4=0.1;
Ro=1000;
% G function and C function
% G*V+jwC*V=F
% G=[R3 -1 0 0;R3 0 -1 0;alpha*Ro/(R4+Ro) 0 0 -1;1 0 (1/R1+1/R2) 0];
G=[R3 -1 0 0 0;R3 0 -1 0 0;alpha*Ro/(R4+Ro) 0 0 -1 0;1 0 (1/R1+1/R2) 0 0;0 0 0 0 1];
% C=[0 0 0 0;L 0 0 0;0 0 0 0;0 0 c 0];
C=[0 0 0 0 0;L 0 0 0 0;0 0 0 0 0;0 0 c 0 -c/R1;0 0 0 0 0];
dt=1/1000; % time interval
A=C/dt+G; % prepare for time stepping

% a: step input
t=0; % initial time
Vina=zeros(1,1000); % initialize size for input
Voa=zeros(1,1000); % initialize size for output
for j=1:1:1000
    if t<0.03
        Vina(j)=0;
    else
        Vina(j)=1;
    end
    FDC=[0;0;0;Vina(j)/R1;Vina(j)]; % corresponding F
    if j==1
        DC=zeros(5,1); % V(0)
    else
        DC=A\(C*oldV/dt+FDC); % A*V(n)=C*V(n-1)/dt+F
    end
 
    Voa(j)=DC(4);
    oldV=DC; % update V(n-1)
    t=t+dt;
end
figure(1)
t=linspace(0,1,1000);
subplot(1,2,1),plot(t,Vina);
title('Vin(a) vs. t');
xlabel('t');
ylabel('Vin(a)');
grid on
subplot(1,2,2),plot(t,Voa);
title('Vo(a) vs. t');
xlabel('t');
ylabel('Vo(a)');
grid on 

n1=2^nextpow2(1000);
m1=fftshift(fft(Vina,n1+1));
m2=fftshift(fft(Voa,n1+1));
f1=1/0.03*((-n1/2):(n1/2))/n1;

figure(2)
subplot(1,2,1),semilogy(f1,abs(m1/n1));
title('Vin(a) in frequency domain');
xlabel('freq');
grid on
subplot(1,2,2),semilogy(f1,abs(m2/n1));
title('Vo(a) in frequency domain');
xlabel('freq');
grid on 

% b
t=0;
f=1/0.03;
w=2*pi*f;
Vinb=zeros(1,1000);
Vob=zeros(1,1000);
for j=1:1:1000
    Vinb(j)=sin(2*pi*f*t);
    F=[0;0;0;Vinb(j)/R1+c*w*1i*Vinb(j);Vinb(j)];
    if j==1
        V=zeros(5,1); 
    else
        V=A\(C*oldVb/dt+F);
    end
    Vob(j)=abs(V(4));
    oldVb=V;
    t=t+dt;
end
figure(3)
t=linspace(0,1,1000);
subplot(1,2,1),plot(t,sin(2*pi*f*t),'g');
title('Vin(b) vs. t');
xlabel('t');
ylabel('Vin(b)');
grid on
subplot(1,2,2),plot(t,Vob,'b');
title('Vo(b) vs. t');
xlabel('t');
ylabel('Vo(b)');
grid on 

n2=2^nextpow2(1000);
dim=2;
m3=fftshift(fft(Vinb,n2-1,dim));
m4=fftshift(fft(Vob,n2-1,dim));
f2=(f/n2-f/2):f/n2:(f/2-f/n2);
figure(4)
subplot(1,2,1),semilogy(f2,abs(m3/1000));
title('Vin(b) in frequency domain');
xlabel('freq');
grid on
subplot(1,2,2),semilogy(f2,abs(m4/1000));
title('Vo(b) in frequency domain');
xlabel('freq');
grid on 

% c
t=0;
Vinc=zeros(1,1000);
Voc=zeros(1,1000);
Vinc(1:30)=normpdf(0.001:0.001:0.03,0.015,0.03);
[index,m]=max(Vinc(1:30));
M=m/0.06;
Vinc(1:30)=Vinc(1:30)/M;
for k=91:1000
    z=mod(k,90);
    if z==0
        Vinc(k)=0;
    else
        Vinc(k)=Vinc(z);
    end
end
for j=1:1:1000
    FGa=[0;0;0;Vinc(j)/R1+c*w*1i*Vinc(j);Vinc(j)];
    if j==1
        VDC=zeros(5,1); 
    else
        VDC=A\(C*oldVc/dt+FGa);
    end
    Voc(j)=abs(VDC(4));
    oldVc=VDC;
    t=t+dt;
end
figure(5)
t=linspace(0,1,1000);
subplot(1,2,1),plot(t,Vinc);
title('Vin(c) vs. t');
xlabel('t');
ylabel('Vin(c)');
grid on
subplot(1,2,2),plot(t,Voc);
title('Vo(c) vs. t');
xlabel('t');
ylabel('Vo(c)');
grid on 

n3=2^nextpow2(1000);
m5=fftshift(fft(Vinc,n3+1));
m6=fftshift(fft(Voc,n3+1));
f3=1/0.03*(-n3/2:n3/2)/n3;
figure(6)
subplot(1,2,1),semilogy(f3,abs(m5/n3));
title('Vin(c) in frequency domain');
xlabel('freq');
grid on
subplot(1,2,2),semilogy(f3,abs(m6/n3));
title('Vo(c) in frequency domain');
xlabel('freq');
grid on 

