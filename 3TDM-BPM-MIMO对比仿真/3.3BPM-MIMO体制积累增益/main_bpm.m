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
V_Max = radar_parameter.lambda/4/Tp;        % 最大不模糊速度
V_res = V_Max*2/slow_num;                   % 速度分辨率
% 这里仅考虑均匀排布天线阵列
Tx_num = 2;                                     % 发射天线数
Rx_num = 4;                                     % 接收天线数
dRx = lambda/2;
Rx_data = zeros(Rx_num,sample_num,slow_num);
% 发射天线二进制相位编码，这里仅可用于2发形式
bpm_value_1 = kron(ones(sample_num,1),exp(1i*zeros(1,slow_num)));      % 扩展以便和矩阵相乘
bpm_value_2 = kron(ones(sample_num,1),exp(1i*pi*(0:slow_num-1)));
% 接收通道数据
for Rx_id = 1:Rx_num
    ADC_data = zeros(sample_num,slow_num);
    Tx_data_1 = ADC_data;
    Tx_data_2 = ADC_data;
    for Tx_id = 1:Tx_num
        temp_data = zeros(sample_num,slow_num);
        % 叠加所有目标
        for temp_tar = 1:size(tar_inf,1)
            % slow_time = (Tx_id - 1 + Tx_num*(0:slow_num-1))*Tp;                   % TDM模式下时序
            slow_time = (0:slow_num-1)*Tp;                                          % 这里和TDM不同
            tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;          % 生成波形，距离R=R0+v*dt
            fc_array = radar_parameter.fc*ones(sample_num,1);
            fast_time = ((0:sample_num-1)/fs)';
            delay = 2*tar_pos/c;
            phi = ((Tx_id-1)*Rx_num+Rx_id)*dRx/lambda * sind(tar_inf(temp_tar,3));  % 收发天线相位差，仿真DOA
            temp_data = temp_data + exp(1i*2*pi*(phi + k*fast_time*delay + fc_array*delay  - 1/2*k*(delay).^2));
        end
        % 叠加二进制相位编码
        if Tx_id == 1
            Tx_data_1 = temp_data .* bpm_value_1;
        else
            Tx_data_2 = temp_data .* bpm_value_2;
        end
    end
    ADC_data = Tx_data_1 + Tx_data_2;
    ADC_data = awgn(ADC_data,0);        % 产生噪声,0dB
    Rx_data(Rx_id,:,:) = ADC_data;
end
% 显示解调前的数据，包括时域和二维FFT数据
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
[XX,YY] = meshgrid(V_label,R_label);
% 选取第1通道显示
Rx_id = 1;
temp_data = squeeze(Rx_data(Rx_id,:,:));
FFT_1D = fft(temp_data,sample_num,1);
FFT_2D = abs((fftshift(fft(FFT_1D,[],2),2)));
% 各通道解调，解调后chirp周期翻倍，最大不模糊速度下降一半
slow_num = slow_num/2;
V_Max = radar_parameter.lambda/4/(2*Tp);        % 最大不模糊速度下降一半
V_res = V_Max*2/slow_num;                   % 速度分辨率不变
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
[XX,YY] = meshgrid(V_label,R_label);
Channel_data = zeros(Tx_num*Rx_num,sample_num,slow_num);
for Rx_id = 1:Rx_num
    temp_data = squeeze(Rx_data(Rx_id,:,:));
    decode_data_1 = zeros(sample_num, slow_num);
    decode_data_2 = zeros(sample_num, slow_num);
    for i = 1:slow_num
        Sa = temp_data(:, 2*i-1);
        Sb = temp_data(:, 2*i);
        S1 = (Sa + Sb)/2;
        S2 = (Sa - Sb)/2;
        % 调试用，单各chirp的信号
%         figure(1);
%         subplot(221);plot(abs(Sa)); title('Sa信号幅值');grid on;
%         subplot(222);plot(abs(Sb)); title('Sb信号幅值');grid on;
%         subplot(223);plot(abs(S1)); title('S1信号幅值');grid on;
%         subplot(224);plot(abs(S2)); title('S2信号幅值');grid on;
%         figure(2);
%         FFT_1D = abs(fft(Sa,sample_num,1));subplot(221);plot(db(FFT_1D)); title('Sa信号幅值');grid on;
%         FFT_1D = abs(fft(Sb,sample_num,1));subplot(222);plot(db(FFT_1D)); title('Sb信号幅值');grid on; 
%         FFT_1D = abs(fft(S1,sample_num,1));subplot(223);plot(db(FFT_1D)); title('S1信号幅值');grid on; 
%         FFT_1D = abs(fft(S2,sample_num,1));subplot(224);plot(db(FFT_1D)); title('S2信号幅值');grid on;
        decode_data_1(:, i) = S1;
        decode_data_2(:, i) = S2;
    end
    Channel_data(Rx_id,:,:) = decode_data_1;
    Channel_data(Rx_id + Rx_num,:,:) = decode_data_2;
    % 显示各个通道的RD图
    % 一维FFT和二维FFT
    FFT_1D = fft(decode_data_1,sample_num,1);
    FFT_2D_1 = (fftshift(fft(FFT_1D,[],2),2));
    Channel_data(Rx_id,:,:) = FFT_2D_1;
    FFT_1D = fft(decode_data_2,sample_num,1);
    FFT_2D_2 = fftshift(fft(FFT_1D,[],2),2);
    Channel_data(Rx_id + Rx_num,:,:) = FFT_2D_2;
%     figure;mesh(XX,YY,db(abs(FFT_2D_1)));colorbar;title('单通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
%     figure;subplot(121);
%     imagesc(V_label,R_label,db(abs(FFT_2D_1)));colorbar;title('单通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
%     subplot(122);
%     imagesc(V_label,R_label,db(abs(FFT_2D_2)));colorbar;title('单通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
end


% 显示各个通道的RD图
Acc_channel_data = zeros(sample_num,slow_num);
% figure;
for i = 1:Tx_num*Rx_num
%     subplot(2,4,i);
    temp_data = abs(squeeze(Channel_data(i,:,:)));
    Acc_channel_data = Acc_channel_data + temp_data;
end

mesh(V_label,R_label,db(Acc_channel_data));title('积累后二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
max_value = max(max(db(Acc_channel_data)));
noise_value = db(mean(mean(Acc_channel_data)));
snr = max_value - noise_value;

