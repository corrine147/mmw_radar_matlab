% 确定性最大似然估计仿真
clc; clear; close all;

%% 参数设置
%%% 工作频率
lambda = 1;         % 波长
k = 2*pi/lambda;    % 波数
%%% 阵列参数
N = 10;                 % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(N-1)*d)';     % 阵元坐标分布
n_list = (0:N-1)';
%%% 信号源参数
% target_theta = [-30, -10.6, 20.4, 60.2]';     % 多个来波方向
target_theta = [-30, 35]';     % 多个来波方向
M = length(target_theta);           % 信号源数目
%%% 仿真参数
SNR = 10;     % 信噪比(dB)
L = 500;     % 采样点数

%% 阵列接收信号仿真模拟
A = exp(1j*k*z*sind(target_theta'));          % 流型矩阵
S = randn(M, L);         % 输入信号
X = A*S;                        % 阵列接收信号
X1 = awgn(X, SNR, 'measured');      % 加载高斯白噪声后的阵列接收信号
% 阵列接收信号的协方差矩阵
R = X1*X1'/L;    

%% 确定初始值
theta_list = (-90:0.1:90)'*pi/180;
for idx = 1:length(theta_list)
    %%% 计算DML谱函数
    psi0 = k*d*sin(theta_list(idx));        % electrical angle
    A = exp(1j*n_list*psi0');
    P_A = A*inv(A'*A)*A';                       % projection matrix wrt V
    f(idx) = trace(P_A*R);
end
%%% 谱函数归一化
f = real(f);
f = (f-min(f))/(max(f)-min(f));
[~, ii] = max(f);       % 峰值
theta0 = theta_list(ii);
psi = k*d*sin(theta0);      % 将峰值作为牛顿迭代初始值
%%% 获取谱峰值
[f_peaks, f_peaks_idx] = findpeaks(f);     % 提取峰值
[f_peaks, I] = sort(f_peaks, 'descend');    % 峰值降序排序
f_peaks_idx = f_peaks_idx(I);
f_peaks = f_peaks(1:M);             % 提取前M个
f_peaks_idx = f_peaks_idx(1:M);
theta_est = theta_list(f_peaks_idx)*180/pi;   % 估计方向
disp('信号源估计方向为：');
disp(theta_est);
%%% 绘制 DML 谱函数
figure;
plot(theta_list*180/pi, db(f));hold on;
xlabel('\theta (deg)');ylabel('DML空间谱');title('确定性最大似然估计角度');
plot([target_theta,target_theta]',ylim,'m-.');hold on;
plot(theta_est, db(f_peaks), 'r.', 'MarkerSize', 15);hold on;
for idx = 1:M
    text(theta_est(idx)+3, f_peaks(idx), sprintf('%0.1f°', theta_est(idx)));
end
grid on;ylim([-90,10]);