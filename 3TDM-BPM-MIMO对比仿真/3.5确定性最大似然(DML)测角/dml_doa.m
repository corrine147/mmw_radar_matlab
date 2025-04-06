% DML用于角度估计，代码有点问题
clear;close all;clc;
 
nT = 8; % 发射阵元
nR = 8; % 接收阵元
L = 181;
M = 101;
theta = linspace(-90,90,L).'*pi/180;
freq = linspace(-0.5,0.5,M).';
 
% 发送信号模型
N = 100;
X = (sign(randn(nT,N)) + 1j*sign(randn(nT,N)))/sqrt(2);
varn = 1;
 
% 目标2参数
fd = -0.1;
thetat = -55*pi/180;
d = exp(1j*2*pi*fd*(0:N-1)).';
aT = exp(1j*pi*(0:nT-1)*sin(thetat)).';
aR = exp(1j*pi*(0:nR-1)*sin(thetat)).';
Y = aR*aT.'*X*diag(d);
 
% 目标2参数
fd = 0.3;
thetat = -55*pi/180;
d = exp(1j*2*pi*fd*(0:N-1)).';
aT = exp(1j*pi*(0:nT-1)*sin(thetat)).';
aR = exp(1j*pi*(0:nR-1)*sin(thetat)).';
Y = Y + aR*aT.'*X*diag(d);
 
% 高斯噪声
Y = Y + sqrt(varn)*(randn(nR,N) + 1j*randn(nR,N));
y = vec(Y);
 
% Detection
P = zeros(L,M); % L为遍历的角度索引值，M为遍历的多普勒索引值
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