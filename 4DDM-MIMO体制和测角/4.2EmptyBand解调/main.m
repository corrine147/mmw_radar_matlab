% DDM体制FMCW雷达目标回波、FFT、FFT_1D、FFT_2D仿真
clc;clear;close all;
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
Tx_num = 4;                                     % 发射天线数
EmptyBand_num = 2;                              % 空带数量
Total_num = Tx_num + EmptyBand_num;             % 总的发射数量为实际发射天线加上空带
Rx_num = 4;                                     % 接收天线数
dRx = lambda/2;
Rx_data = zeros(Rx_num,sample_num,slow_num);
% 发射天线相位编码，理论上可用于任意根天线
ddm_phase_value_1 = kron(ones(sample_num,1),exp(1i*zeros(1,slow_num)));      % 扩展以便和矩阵相乘
ddm_phase_value_2 = kron(ones(sample_num,1),exp(1i*2*pi*(0:slow_num-1)/Total_num));
ddm_phase_value_3 = kron(ones(sample_num,1),exp(1i*4*pi*(0:slow_num-1)/Total_num));
ddm_phase_value_4 = kron(ones(sample_num,1),exp(1i*6*pi*(0:slow_num-1)/Total_num));
% 目标参数
tar_inf = [55 6+2*V_Max/Total_num 0];          % 仿真目标距离,速度,角度
% tar_inf = [15 1 10;
%            18 -5 0;
%            70.8 20 0;
%            90.3 -15 20;
%            105 10 -30;];          % 仿真目标距离,速度,角度
% 接收通道数据
for Rx_id = 1:Rx_num
    % 每根发射天线的数据
    ADC_data = zeros(sample_num,slow_num);
    Tx_data_1 = ADC_data;
    Tx_data_2 = ADC_data;
    Tx_data_3 = ADC_data;
    Tx_data_4 = ADC_data;
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
        % 叠加慢时间相位编码
        switch Tx_id
            case 1
                Tx_data_1 = temp_data .* ddm_phase_value_1;
            case 2
                Tx_data_2 = temp_data .* ddm_phase_value_2;
            case 3
                Tx_data_3 = temp_data .* ddm_phase_value_3;
            case 4
                Tx_data_4 = temp_data .* ddm_phase_value_4;
            otherwise
                disp('default');
        end
    end
%     ADC_data = Tx_data_3;             % 单根天线发射用于测试
    ADC_data = Tx_data_1 + Tx_data_2 + Tx_data_3 + Tx_data_4; % 所有天线发射
    ADC_data = awgn(ADC_data,0);        % 产生噪声,0dB
    Rx_data(Rx_id,:,:) = ADC_data;
end
% 显示解调前的数据，包括时域和二维FFT数据
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
[XX,YY] = meshgrid(V_label,R_label);
% 选取第1接收通道显示
Rx_id = 1;
temp_data = squeeze(Rx_data(Rx_id,:,:));
FFT_1D = fft(temp_data,sample_num,1);
FFT_2D = abs((fftshift(fft(FFT_1D,[],2),2)));
figure;mesh(XX,YY,db(abs(FFT_1D)));colorbar;title('第1通道一维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;mesh(XX,YY,db(FFT_2D));colorbar;title('第1通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
% figure;imagesc(V_label,R_label,db(FFT_2D));colorbar;title('第1单通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
Image_fft = db(FFT_2D);             % 将最大不模糊速度根据发射天线数划分
V_label_len = length(V_label);      % 速度范围的chirp数，根据划分定位各通道分界线
for i = 1:Total_num-1
    Line = round(V_label_len/Total_num)*i;
    Image_fft(:,Line) = 0;
end
figure;imagesc(V_label,R_label,Image_fft);colorbar;title('第1单通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');




