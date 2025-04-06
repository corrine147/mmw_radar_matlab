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
V_Max = radar_parameter.lambda/4/Tp;        % ���ģ���ٶ�
V_res = V_Max*2/slow_num;                   % �ٶȷֱ���
% ��������Ǿ����Ų���������
Tx_num = 2;                                     % ����������
Rx_num = 4;                                     % ����������
dRx = lambda/2;
Rx_data = zeros(Rx_num,sample_num,slow_num);
% �������߶�������λ���룬�����������2����ʽ
bpm_value_1 = kron(ones(sample_num,1),exp(1i*zeros(1,slow_num)));      % ��չ�Ա�;������
bpm_value_2 = kron(ones(sample_num,1),exp(1i*pi*(0:slow_num-1)));
% ����ͨ������
for Rx_id = 1:Rx_num
    ADC_data = zeros(sample_num,slow_num);
    Tx_data_1 = ADC_data;
    Tx_data_2 = ADC_data;
    for Tx_id = 1:Tx_num
        temp_data = zeros(sample_num,slow_num);
        % ��������Ŀ��
        for temp_tar = 1:size(tar_inf,1)
            % slow_time = (Tx_id - 1 + Tx_num*(0:slow_num-1))*Tp;                   % TDMģʽ��ʱ��
            slow_time = (0:slow_num-1)*Tp;                                          % �����TDM��ͬ
            tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;          % ���ɲ��Σ�����R=R0+v*dt
            fc_array = radar_parameter.fc*ones(sample_num,1);
            fast_time = ((0:sample_num-1)/fs)';
            delay = 2*tar_pos/c;
            phi = ((Tx_id-1)*Rx_num+Rx_id)*dRx/lambda * sind(tar_inf(temp_tar,3));  % �շ�������λ�����DOA
            temp_data = temp_data + exp(1i*2*pi*(phi + k*fast_time*delay + fc_array*delay  - 1/2*k*(delay).^2));
        end
        % ���Ӷ�������λ����
        if Tx_id == 1
            Tx_data_1 = temp_data .* bpm_value_1;
        else
            Tx_data_2 = temp_data .* bpm_value_2;
        end
    end
    ADC_data = Tx_data_1 + Tx_data_2;
    ADC_data = awgn(ADC_data,0);        % ��������,0dB
    Rx_data(Rx_id,:,:) = ADC_data;
end
% ��ʾ���ǰ�����ݣ�����ʱ��Ͷ�άFFT����
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
[XX,YY] = meshgrid(V_label,R_label);
% ѡȡ��1ͨ����ʾ
Rx_id = 1;
temp_data = squeeze(Rx_data(Rx_id,:,:));
FFT_1D = fft(temp_data,sample_num,1);
FFT_2D = abs((fftshift(fft(FFT_1D,[],2),2)));
% ��ͨ������������chirp���ڷ��������ģ���ٶ��½�һ��
slow_num = slow_num/2;
V_Max = radar_parameter.lambda/4/(2*Tp);        % ���ģ���ٶ��½�һ��
V_res = V_Max*2/slow_num;                   % �ٶȷֱ��ʲ���
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
        % �����ã�����chirp���ź�
%         figure(1);
%         subplot(221);plot(abs(Sa)); title('Sa�źŷ�ֵ');grid on;
%         subplot(222);plot(abs(Sb)); title('Sb�źŷ�ֵ');grid on;
%         subplot(223);plot(abs(S1)); title('S1�źŷ�ֵ');grid on;
%         subplot(224);plot(abs(S2)); title('S2�źŷ�ֵ');grid on;
%         figure(2);
%         FFT_1D = abs(fft(Sa,sample_num,1));subplot(221);plot(db(FFT_1D)); title('Sa�źŷ�ֵ');grid on;
%         FFT_1D = abs(fft(Sb,sample_num,1));subplot(222);plot(db(FFT_1D)); title('Sb�źŷ�ֵ');grid on; 
%         FFT_1D = abs(fft(S1,sample_num,1));subplot(223);plot(db(FFT_1D)); title('S1�źŷ�ֵ');grid on; 
%         FFT_1D = abs(fft(S2,sample_num,1));subplot(224);plot(db(FFT_1D)); title('S2�źŷ�ֵ');grid on;
        decode_data_1(:, i) = S1;
        decode_data_2(:, i) = S2;
    end
    Channel_data(Rx_id,:,:) = decode_data_1;
    Channel_data(Rx_id + Rx_num,:,:) = decode_data_2;
    % ��ʾ����ͨ����RDͼ
    % һάFFT�Ͷ�άFFT
    FFT_1D = fft(decode_data_1,sample_num,1);
    FFT_2D_1 = (fftshift(fft(FFT_1D,[],2),2));
    Channel_data(Rx_id,:,:) = FFT_2D_1;
    FFT_1D = fft(decode_data_2,sample_num,1);
    FFT_2D_2 = fftshift(fft(FFT_1D,[],2),2);
    Channel_data(Rx_id + Rx_num,:,:) = FFT_2D_2;
%     figure;mesh(XX,YY,db(abs(FFT_2D_1)));colorbar;title('��ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
%     figure;subplot(121);
%     imagesc(V_label,R_label,db(abs(FFT_2D_1)));colorbar;title('��ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
%     subplot(122);
%     imagesc(V_label,R_label,db(abs(FFT_2D_2)));colorbar;title('��ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
end


% ��ʾ����ͨ����RDͼ
Acc_channel_data = zeros(sample_num,slow_num);
% figure;
for i = 1:Tx_num*Rx_num
%     subplot(2,4,i);
    temp_data = abs(squeeze(Channel_data(i,:,:)));
    Acc_channel_data = Acc_channel_data + temp_data;
end

mesh(V_label,R_label,db(Acc_channel_data));title('���ۺ��άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
max_value = max(max(db(Acc_channel_data)));
noise_value = db(mean(mean(Acc_channel_data)));
snr = max_value - noise_value;

