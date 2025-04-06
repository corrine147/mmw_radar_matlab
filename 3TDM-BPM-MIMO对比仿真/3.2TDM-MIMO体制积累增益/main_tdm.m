% FMCW目标回波、FFT、FFT_1D、FFT_2D、CFAR仿真
clc;clear;close all;
% 目标参数
tar_inf = [50 0 0;];          % 仿真目标距离,速度,角度
% 雷达参数设置
radar_parameter = signal_para_set();
c = radar_parameter.c;
k = radar_parameter.K;
Tp = radar_parameter.Tp;
fs = radar_parameter.fs;
lambda = radar_parameter.lambda;
sample_num = radar_parameter.sample_num;        	% 单chirp采样点数
slow_num = radar_parameter.slow_num;                % chirp数量

R_res = c/2/radar_parameter.real_B;             % 距离分辨率
R_Max = R_res*sample_num;                       % 最远距离
V_Max = radar_parameter.lambda/4/(2*Tp);        % 最大不模糊速度,相比于单发模式降低了一半
V_res = V_Max*2/slow_num;                       % 速度分辨率，相比于单发模式也降低了一半
% 这里仅考虑均匀排布天线阵列
Tx_num = 2;                                     % 发射天线数
Rx_num = 4;                                     % 接收天线数
dRx = lambda/2;
Channel_data = zeros(Tx_num*Rx_num,sample_num,slow_num);
for Tx_id = 1:Tx_num
    for Rx_id = 1:Rx_num
        ADC_data = zeros(sample_num,slow_num);
        for temp_tar = 1:size(tar_inf,1)
            slow_time = (Tx_id - 1 + 2*(0:slow_num-1))*Tp;
            tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;        % 生成波形，距离R=R0+v*dt
            % temp_data = signal_generate(1, 2*tar_pos, radar_parameter);
            fc_array = radar_parameter.fc*ones(sample_num,1);
            fast_time = ((0:sample_num-1)/fs)';
            delay = 2*tar_pos/c;
            phi = ((Tx_id-1)*Rx_num+Rx_id)*dRx/lambda * sind(tar_inf(temp_tar,3));
            temp_data = exp(1i*2*pi*(phi + k*fast_time*delay + fc_array*delay  - 1/2*k*(delay).^2));
            ADC_data = ADC_data + temp_data;
        end
        ADC_data = awgn(ADC_data,0);        % 产生噪声,0dB
        % 1D FFT
        FFT_1D = fft(ADC_data,sample_num,1);
        % 2D FFT
        FFT_2D = fftshift(fft(FFT_1D,[],2),2);
        Channel_data((Tx_id-1)*Rx_num+Rx_id,:,:) = FFT_2D;
    end
end

% 显示各个通道的RD图
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
Acc_channel_data = zeros(sample_num,slow_num);
% figure;
for i = 1:Tx_num*Rx_num
%     subplot(2,4,i);
    temp_data = abs(squeeze(Channel_data(i,:,:)));
    Acc_channel_data = Acc_channel_data + temp_data;
end
mesh(V_label,R_label,db(Acc_channel_data));title('积累后二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');

% 单通道信噪比



% 多通道积累后信噪比
max_value = max(max(db(Acc_channel_data)));
noise_value = db(mean(mean(Acc_channel_data)));
snr = max_value - noise_value;








