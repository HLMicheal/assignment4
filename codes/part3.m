% ELEC4700 Assignment 4
% Part 3
% By Huanyu Liu
% 100986552
clear
clc
% 
CN=[0.00001,0.001,0.00000005,0.00001,0.00001];
DT=[1/1000,1/1000,1/1000,1/5000,1/2000];
    R1=1;
    c=0.25;
    R2=2;
    L=0.2;
    R3=10;
    alpha=100;
    R4=0.1;
    Ro=1000;
for k=1:5
    
    cn=CN(k);
    G=[R3 0 -1 0 0 0;R3 0 0 -1 0 0;alpha*Ro/(R4+Ro) 0 0 0 -1 0;0 1 0 (1/R1+1/R2) 0 0;0 0 0 0 0 1];
    C=[0 0 0 0 0 0;0 L 0 0 0 0;0 0 0 0 0 0;0 0 0 c 0 -c/R1;0 0 0 0 0 0];

%     G=[R3 0 -1 0 0 0;alpha*Ro/(R4+Ro) 0 0 0 -1 0;0 1 0 1/R1+1/R2 0 0;0 0 -1 1 0 0;-1 1 0 0 0 0];
%     C=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 c 0 -c/R1;0 -L 0 0 0 0;0 0 cn 0 0 0];
    % % % V=[I3;Iall;V3;V2;Vo;V1];

    t=0;
    dt=DT(k);
    f=1/0.03;
    w=2*pi*f;
    Vinc=zeros(1,1000);
    Voc=zeros(1,1000);
    Vinc(1:30)=normpdf(0.001:0.001:0.03,0.015,0.03);
    [index,m]=max(Vinc(1:30));
    M=m/0.06;
    Vinc(1:30)=Vinc(1:30)/M;
    for m=91:1000
        z=mod(m,90);
        if z==0
            Vinc(m)=0;
        else
            Vinc(m)=Vinc(z);
        end
    end
    for j=1:1:1000
        In=normrnd(0.001,0.001); %assume std. deviation is 0.001
        FGa=[0;0;Vinc(j)/R1+c*w*1i*Vinc(j);0;In];
        if j==1
            VDC=zeros(6,1); %G\FGa;
        else
            VDC=A\(C*oldV/dt+FGa);
        end
        Voc(j)=abs(VDC(4));
        oldV=VDC;
        A=C/dt+G;
        t=t+dt;
    end
    figure(2*k-1)
    t=linspace(0,1,1000);
    subplot(1,2,1),plot(t,Vinc,'g');
    title(['Vin vs. t (C=',num2str(cn),'dt=',num2str(dt),')']);
    xlabel('t');
    ylabel('Vin(part3)');
    grid on
    subplot(1,2,2),plot(t,Voc,'b');
    title(['Vo vs. t (C=',num2str(cn),'dt=',num2str(dt),')']);
    xlabel('t');
    ylabel('Vo(part3)');
    grid on 

    n=2^nextpow2(1000);
    m5=fftshift(fft(Vinc,n+1));
    m6=fftshift(fft(Voc,n+1));
    f=1/0.03*((-n/2):(n/2))/n;
    figure(2*k)
    subplot(1,2,1),plot(f,abs(m5/n),'r');
    title(['Vin in frequency domain (C=',num2str(cn),'dt=',num2str(dt),')']);
    xlabel('freq');
    grid on
    subplot(1,2,2),plot(f,abs(m6/n),'b');
    title(['Vo in frequency domain (C=',num2str(cn),'dt=',num2str(dt),')']);
    xlabel('freq');
    grid on 
end