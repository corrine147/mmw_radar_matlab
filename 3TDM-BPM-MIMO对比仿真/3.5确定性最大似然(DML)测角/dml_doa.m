% DML���ڽǶȹ��ƣ������е�����
clear;close all;clc;
 
nT = 8; % ������Ԫ
nR = 8; % ������Ԫ
L = 181;
M = 101;
theta = linspace(-90,90,L).'*pi/180;
freq = linspace(-0.5,0.5,M).';
 
% �����ź�ģ��
N = 100;
X = (sign(randn(nT,N)) + 1j*sign(randn(nT,N)))/sqrt(2);
varn = 1;
 
% Ŀ��2����
fd = -0.1;
thetat = -55*pi/180;
d = exp(1j*2*pi*fd*(0:N-1)).';
aT = exp(1j*pi*(0:nT-1)*sin(thetat)).';
aR = exp(1j*pi*(0:nR-1)*sin(thetat)).';
Y = aR*aT.'*X*diag(d);
 
% Ŀ��2����
fd = 0.3;
thetat = -55*pi/180;
d = exp(1j*2*pi*fd*(0:N-1)).';
aT = exp(1j*pi*(0:nT-1)*sin(thetat)).';
aR = exp(1j*pi*(0:nR-1)*sin(thetat)).';
Y = Y + aR*aT.'*X*diag(d);
 
% ��˹����
Y = Y + sqrt(varn)*(randn(nR,N) + 1j*randn(nR,N));
y = vec(Y);
 
% Detection
P = zeros(L,M); % LΪ�����ĽǶ�����ֵ��MΪ�����Ķ���������ֵ
for i = 1:L
    aT = exp(1j*pi*(0:nT-1)*sin(theta(i))).';
    aR = exp(1j*pi*(0:nR-1)*sin(theta(i))).';
    A = aR*aT.';
    for k = 1:M
        dk = exp(1j*2*pi*freq(k)*(0:N-1)).';
        vk = vec(A*X*diag(dk));
        P(i,k) = abs(y.'*vk + vk.'*y) / (2*N*nT*nR);
    end
end
 
figure;
mesh(freq,theta*180/pi,P)
xlabel('Frequency')
ylabel('DOA (degree)')
zlabel('Post-Correlation Power')
title('Ntx = 8')