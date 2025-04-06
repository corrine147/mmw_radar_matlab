% 幅相不一致对频谱峰值的影响
clc;clear;close all;
N = 16;                 % 通道数
L = 1024;               % 时域信号序列长度      
fs = 6000;              % 采样频率
ts = 1/fs;              % 采样时间间隔
t = ts:ts:L*ts;         % 时间序列
f = 1000;               % 信号频率
phase_error = 6;        % 3,6,9
phase = phase_error * randn(1,N)*pi/180;        % 每个通道的初始相位，单位弧度
A_error = 0.15;         % 0.05,0.10,0.15
A = (1-A_error) + A_error * randn(1,N);         % 每个通道的初始幅度
times = 1500;             % 循环次数
snr_loss = zeros(1,times);
for k = 1:times
    s_clean = zeros(N, L);  % 无噪声信号       
    s_noise = zeros(N, L);  % 有噪声信号   
    % 创建仿真信号
    for i = 1:N
        s_clean(i,:) = cos(2*pi*f*t) + randn(1,L);
        s_noise(i,:) = A(i)*cos(2*pi*f*t + phase(i)) + randn(1,L);
    end
    % 画出时域信号

    Nf = 8192;              % 计算FFT的点数
    df = fs/Nf;              % 频率分辨率
    s_clean_fft = zeros(N,Nf);     % s_clean的FFT结果 
    s_noise_fft = zeros(N,Nf);     % s_noise的FFT结果 
    for i = 1:N
        s_clean_fft(i,:) = abs(fft(s_clean(i,:),Nf));
        s_noise_fft(i,:) = abs(fft(s_noise(i,:),Nf));
    end
    sum_clean_fft = sum(s_clean_fft);
    sum_noise_fft = sum(s_noise_fft);
    [~,clean_snr] = get_max_value(sum_clean_fft, Nf, df);
    [~,noise_snr] = get_max_value(sum_noise_fft, Nf, df);
    snr_loss(k) = clean_snr - noise_snr;
end
snr_loss_mean = mean(snr_loss);
display(snr_loss_mean);

% figure;
% plot(1:Nf/2, db(sum_clean_fft(1:Nf/2)));hold on;
% plot(1:Nf/2, db(sum_noise_fft(1:Nf/2)));hold on;
% grid on;
% legend('幅相一致','幅相不一致');
% title('幅相一致和不一致的对比');

% 通过最大幅值响应计算频率和信噪比
function [f_cal, snr] = get_max_value(sum_fft, Nf, df)
    [fft_max, index] = max(sum_fft,[],2);
    f_cal = index*df;
    not_noise_num = Nf/8;       % 峰值两侧非底噪点数
    max_value = db(fft_max);    
    noise_value = db((sum(sum_fft(1:index - not_noise_num)) + sum(sum_fft(index + not_noise_num:end)))/(3*Nf/4));
    snr = max_value - noise_value;
end







