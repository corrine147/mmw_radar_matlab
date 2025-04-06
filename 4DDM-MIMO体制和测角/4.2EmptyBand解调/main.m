% DDM����FMCW�״�Ŀ��ز���FFT��FFT_1D��FFT_2D����
clc;clear;close all;
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
Tx_num = 4;                                     % ����������
EmptyBand_num = 2;                              % �մ�����
Total_num = Tx_num + EmptyBand_num;             % �ܵķ�������Ϊʵ�ʷ������߼��Ͽմ�
Rx_num = 4;                                     % ����������
dRx = lambda/2;
Rx_data = zeros(Rx_num,sample_num,slow_num);
% ����������λ���룬�����Ͽ��������������
ddm_phase_value_1 = kron(ones(sample_num,1),exp(1i*zeros(1,slow_num)));      % ��չ�Ա�;������
ddm_phase_value_2 = kron(ones(sample_num,1),exp(1i*2*pi*(0:slow_num-1)/Total_num));
ddm_phase_value_3 = kron(ones(sample_num,1),exp(1i*4*pi*(0:slow_num-1)/Total_num));
ddm_phase_value_4 = kron(ones(sample_num,1),exp(1i*6*pi*(0:slow_num-1)/Total_num));
% Ŀ�����
tar_inf = [55 6+2*V_Max/Total_num 0];          % ����Ŀ�����,�ٶ�,�Ƕ�
% tar_inf = [15 1 10;
%            18 -5 0;
%            70.8 20 0;
%            90.3 -15 20;
%            105 10 -30;];          % ����Ŀ�����,�ٶ�,�Ƕ�
% ����ͨ������
for Rx_id = 1:Rx_num
    % ÿ���������ߵ�����
    ADC_data = zeros(sample_num,slow_num);
    Tx_data_1 = ADC_data;
    Tx_data_2 = ADC_data;
    Tx_data_3 = ADC_data;
    Tx_data_4 = ADC_data;
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
        % ������ʱ����λ����
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
%     ADC_data = Tx_data_3;             % �������߷������ڲ���
    ADC_data = Tx_data_1 + Tx_data_2 + Tx_data_3 + Tx_data_4; % �������߷���
    ADC_data = awgn(ADC_data,0);        % ��������,0dB
    Rx_data(Rx_id,:,:) = ADC_data;
end
% ��ʾ���ǰ�����ݣ�����ʱ��Ͷ�άFFT����
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
[XX,YY] = meshgrid(V_label,R_label);
% ѡȡ��1����ͨ����ʾ
Rx_id = 1;
temp_data = squeeze(Rx_data(Rx_id,:,:));
FFT_1D = fft(temp_data,sample_num,1);
FFT_2D = abs((fftshift(fft(FFT_1D,[],2),2)));
figure;mesh(XX,YY,db(abs(FFT_1D)));colorbar;title('��1ͨ��һάFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;mesh(XX,YY,db(FFT_2D));colorbar;title('��1ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
% figure;imagesc(V_label,R_label,db(FFT_2D));colorbar;title('��1��ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
Image_fft = db(FFT_2D);             % �����ģ���ٶȸ��ݷ�������������
V_label_len = length(V_label);      % �ٶȷ�Χ��chirp�������ݻ��ֶ�λ��ͨ���ֽ���
for i = 1:Total_num-1
    Line = round(V_label_len/Total_num)*i;
    Image_fft(:,Line) = 0;
end
figure;imagesc(V_label,R_label,Image_fft);colorbar;title('��1��ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');




