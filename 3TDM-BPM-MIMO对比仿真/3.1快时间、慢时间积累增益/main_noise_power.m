% FMCWĿ��ز���FFT��FFT_1D��FFT_2D����
clc;clear;close all;
c = 3e8;                        % ����
fc = 76e9;                      % ��ʼƵ��
K  = 30e12;                     % ��Ƶб��
Tp = 20e-6;                     % ��chirpʱ��
B = K*Tp;                       % ʵ�ʴ���
lamda = c/fc;                   % ����
fs = 30e6;                      % ��Ƶ����Ƶ��
ts = 1/fs;                      % ����ʱ����
sample_num = 384;               % ��chirp��������
Ts = ts*sample_num;             % ����ʱ��
slow_num = 512;                 % chirp����
tar_inf = [50 0;];            % ����Ŀ�����,�ٶ�
ADC_data = zeros(sample_num,slow_num);
slow_time = (0:slow_num-1)*Tp;
real_B = Ts*K;
R_res = c/2/real_B;             % ����ֱ���
R_MAX = R_res*sample_num;       % ��Զ����
V_Max = lamda/4/Tp;             % ���ģ���ٶ�
V_res = V_Max*2/slow_num;       % �ٶȷֱ���
% �ṹ���������
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

% ���ɲ��Σ�����R=R0+v*dt
for temp_tar = 1:size(tar_inf,1)
    tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;
    temp_data = signal_generate(1, 2*tar_pos, radar_parameter);
    ADC_data = ADC_data + temp_data;
end

R_label = R_res*(1:sample_num);
rec_snr = -15:1:5;
L = length(rec_snr);
snr = zeros(1,L);
N = 20;                                            % ѭ��N��
for i = 1:L
    tmp_snr = zeros(1,N);
    for j = 1:N
        tmp_adc = awgn(ADC_data,rec_snr(i));        % ��������
        % 1D FFT
        FFT_1D = fft(tmp_adc,sample_num,1);
        % ����һάFFT������ֵ��������
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
xlabel('��������');
ylabel('һάFFT�������-dB');
title('һάFFT��������Ⱥͽ����ź�����ȹ�ϵ');
grid on;


% 2D FFT
% FFT_2D = fftshift(fft(FFT_1D,[],2),2);
% V_lable = (-V_Max:V_res:(V_Max-V_res));
% [XX,YY] = meshgrid(V_lable,R_label);
% % �����άFFT������ֵ��������
% max_value = max(max(db(abs(FFT_2D))));
% noise_value = db(mean(mean(abs(FFT_2D))));
% snr = max_value - noise_value;
% display(max_value);
% display(noise_value);
% display(snr);







