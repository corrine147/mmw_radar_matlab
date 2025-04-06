% FMCW目标回波、FFT、FFT_1D、FFT_2D仿真
clc;clear;close all;
c = 3e8;                        % 光速
fc = 76e9;                      % 起始频率
K  = 30e12;                     % 调频斜率
Tp = 20e-6;                     % 单chirp时长
B = K*Tp;                       % 实际带宽
lamda = c/fc;                   % 波长
fs = 30e6;                      % 中频采样频率
ts = 1/fs;                      % 采样时间间隔
sample_num = 384;               % 单chirp采样点数
Ts = ts*sample_num;             % 采样时间
slow_num = 512;                 % chirp数量
tar_inf = [50 0;];            % 仿真目标距离,速度
ADC_data = zeros(sample_num,slow_num);
slow_time = (0:slow_num-1)*Tp;
real_B = Ts*K;
R_res = c/2/real_B;             % 距离分辨率
R_MAX = R_res*sample_num;       % 最远距离
V_Max = lamda/4/Tp;             % 最大不模糊速度
V_res = V_Max*2/slow_num;       % 速度分辨率
% 结构体参数设置
radar_parameter.c = c;
radar_parameter.B = B;
radar_parameter.Tp = Tp;
radar_parameter.fc = fc;
radar_parameter.K = K;
radar_parameter.fs = fs;
radar_parameter.ts = ts;
radar_parameter.lambda = lamda;
radar_parameter.sample_num = sample_num;
radar_parameter.slow_num = slow_num;
radar_parameter.real_B = real_B;
radar_parameter.R_res = R_res;
radar_parameter.V_Max = V_Max;
radar_parameter.V_res = V_res;

% 生成波形，距离R=R0+v*dt
for temp_tar = 1:size(tar_inf,1)
    tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;
    temp_data = signal_generate(1, 2*tar_pos, radar_parameter);
    ADC_data = ADC_data + temp_data;
end

R_label = R_res*(1:sample_num);
rec_snr = -15:1:5;
L = length(rec_snr);
snr = zeros(1,L);
N = 20;                                            % 循环N次
for i = 1:L
    tmp_snr = zeros(1,N);
    for j = 1:N
        tmp_adc = awgn(ADC_data,rec_snr(i));        % 产生噪声
        % 1D FFT
        FFT_1D = fft(tmp_adc,sample_num,1);
        % 计算一维FFT后最大峰值点的信噪比
        [max_value, max_index] = max(db(abs(FFT_1D(:,1))));
        noise_value = db(mean(abs(FFT_1D(:,1))));
        tmp_snr(j) = max_value - noise_value;
    end
%     if mod(i,2) == 0
%         figure;
%         plot(R_label,db(abs(FFT_1D(:,1))));
%         grid on;
%     end
    snr(i) = mean(tmp_snr);
end

figure;
plot(rec_snr, snr, 'd-');
xlabel('采样点数');
ylabel('一维FFT后信噪比-dB');
title('一维FFT积累信噪比和接收信号信噪比关系');
grid on;


% 2D FFT
% FFT_2D = fftshift(fft(FFT_1D,[],2),2);
% V_lable = (-V_Max:V_res:(V_Max-V_res));
% [XX,YY] = meshgrid(V_lable,R_label);
% % 计算二维FFT后最大峰值点的信噪比
% max_value = max(max(db(abs(FFT_2D))));
% noise_value = db(mean(mean(abs(FFT_2D))));
% snr = max_value - noise_value;
% display(max_value);
% display(noise_value);
% display(snr);







