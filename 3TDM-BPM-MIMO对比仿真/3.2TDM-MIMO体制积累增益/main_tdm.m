% FMCWĿ��ز���FFT��FFT_1D��FFT_2D��CFAR����
clc;clear;close all;
% Ŀ�����
tar_inf = [50 0 0;];          % ����Ŀ�����,�ٶ�,�Ƕ�
% �״��������
radar_parameter = signal_para_set();
c = radar_parameter.c;
k = radar_parameter.K;
Tp = radar_parameter.Tp;
fs = radar_parameter.fs;
lambda = radar_parameter.lambda;
sample_num = radar_parameter.sample_num;        	% ��chirp��������
slow_num = radar_parameter.slow_num;                % chirp����

R_res = c/2/radar_parameter.real_B;             % ����ֱ���
R_Max = R_res*sample_num;                       % ��Զ����
V_Max = radar_parameter.lambda/4/(2*Tp);        % ���ģ���ٶ�,����ڵ���ģʽ������һ��
V_res = V_Max*2/slow_num;                       % �ٶȷֱ��ʣ�����ڵ���ģʽҲ������һ��
% ��������Ǿ����Ų���������
Tx_num = 2;                                     % ����������
Rx_num = 4;                                     % ����������
dRx = lambda/2;
Channel_data = zeros(Tx_num*Rx_num,sample_num,slow_num);
for Tx_id = 1:Tx_num
    for Rx_id = 1:Rx_num
        ADC_data = zeros(sample_num,slow_num);
        for temp_tar = 1:size(tar_inf,1)
            slow_time = (Tx_id - 1 + 2*(0:slow_num-1))*Tp;
            tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;        % ���ɲ��Σ�����R=R0+v*dt
            % temp_data = signal_generate(1, 2*tar_pos, radar_parameter);
            fc_array = radar_parameter.fc*ones(sample_num,1);
            fast_time = ((0:sample_num-1)/fs)';
            delay = 2*tar_pos/c;
            phi = ((Tx_id-1)*Rx_num+Rx_id)*dRx/lambda * sind(tar_inf(temp_tar,3));
            temp_data = exp(1i*2*pi*(phi + k*fast_time*delay + fc_array*delay  - 1/2*k*(delay).^2));
            ADC_data = ADC_data + temp_data;
        end
        ADC_data = awgn(ADC_data,0);        % ��������,0dB
        % 1D FFT
        FFT_1D = fft(ADC_data,sample_num,1);
        % 2D FFT
        FFT_2D = fftshift(fft(FFT_1D,[],2),2);
        Channel_data((Tx_id-1)*Rx_num+Rx_id,:,:) = FFT_2D;
    end
end

% ��ʾ����ͨ����RDͼ
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
Acc_channel_data = zeros(sample_num,slow_num);
% figure;
for i = 1:Tx_num*Rx_num
%     subplot(2,4,i);
    temp_data = abs(squeeze(Channel_data(i,:,:)));
    Acc_channel_data = Acc_channel_data + temp_data;
end
mesh(V_label,R_label,db(Acc_channel_data));title('���ۺ��άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');

% ��ͨ�������



% ��ͨ�����ۺ������
max_value = max(max(db(Acc_channel_data)));
noise_value = db(mean(mean(Acc_channel_data)));
snr = max_value - noise_value;








