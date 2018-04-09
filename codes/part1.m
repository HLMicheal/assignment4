% ELEC4700
% Assignment 4
% Part 1
% By Huanyu Liu
% 100986552
R1=1;
c=0.25;
R2=2;
L=0.2;
R3=10;
alpha=100;
R4=0.1;
Ro=1000;
G=[R3 -1 0 0 0;R3 0 -1 0 0;alpha*Ro/(R4+Ro) 0 0 -1 0;1 0 (1/R1+1/R2) 0 0;0 0 0 0 1];
C=[0 0 0 0 0;L 0 0 0 0;0 0 0 0 0;0 0 c 0 -c/R1;0 0 0 0 0];

figure(1)
title('plot of DC sweep');
for V1=-10:1:10
    FDC=[0;0;0;V1/R1;V1];
    DC=G\FDC;
    plot(V1,DC(4),'r*');
    hold on
    plot(V1,DC(2),'b*');
    hold on
end
xlabel('V1');
ylabel('Vo (red) & V3 (blue)');

figure(2)
title('plots from AC case of gain (Vin=1)');
for w=logspace(-2,4,1000)
    FAC=[0;0;0;1/R1+c*w*1i;1];
    left=G+C*w*1i;
    AC=left\FAC;
    Vo=abs(AC(4));
    subplot(1,2,1),semilogx(w,Vo,'g*');
    title('Vo vs. log10(w)');
    hold on
    grid on
    subplot(1,2,2),semilogx(w,20*log10(Vo),'b*');
    title('(Vo/V1)dB vs. log10(w)');
    hold on
    grid on
end

w=pi;
n=500;
iC=zeros(1,n);
gain=zeros(1,n);
for m=1:n
    iC(m)=normrnd(0.25,0.05);
    nC=[0 0 0 0 0;L 0 0 0 0;0 0 0 0 0;0 0 iC(m) 0 -iC(m)/R1;0 0 0 0 0];
    nFAC=[0;0;0;1+iC(m)*w*1i;1];
    nleft=G+nC*w*1i;
    nAC=nleft\nFAC;
    nVo=abs(nAC(4));
    gain(m)=20*log10(nVo);
end
figure(3)
title('histogram of gain');
histogram(gain);