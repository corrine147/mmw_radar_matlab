% ���಻һ�¶�Ƶ�׷�ֵ��Ӱ��
clc;clear;close all;
N = 16;                 % ͨ����
L = 1024;               % ʱ���ź����г���      
fs = 6000;              % ����Ƶ��
ts = 1/fs;              % ����ʱ����
t = ts:ts:L*ts;         % ʱ������
f = 1000;               % �ź�Ƶ��
phase_error = 6;        % 3,6,9
phase = phase_error * randn(1,N)*pi/180;        % ÿ��ͨ���ĳ�ʼ��λ����λ����
A_error = 0.15;         % 0.05,0.10,0.15
A = (1-A_error) + A_error * randn(1,N);         % ÿ��ͨ���ĳ�ʼ����
times = 1500;             % ѭ������
snr_loss = zeros(1,times);
for k = 1:times
    s_clean = zeros(N, L);  % �������ź�       
    s_noise = zeros(N, L);  % �������ź�   
    % ���������ź�
    for i = 1:N
        s_clean(i,:) = cos(2*pi*f*t) + randn(1,L);
        s_noise(i,:) = A(i)*cos(2*pi*f*t + phase(i)) + randn(1,L);
    end
    % ����ʱ���ź�

    Nf = 8192;              % ����FFT�ĵ���
    df = fs/Nf;              % Ƶ�ʷֱ���
    s_clean_fft = zeros(N,Nf);     % s_clean��FFT��� 
    s_noise_fft = zeros(N,Nf);     % s_noise��FFT��� 
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
% legend('����һ��','���಻һ��');
% title('����һ�ºͲ�һ�µĶԱ�');

% ͨ������ֵ��Ӧ����Ƶ�ʺ������
function [f_cal, snr] = get_max_value(sum_fft, Nf, df)
    [fft_max, index] = max(sum_fft,[],2);
    f_cal = index*df;
    not_noise_num = Nf/8;       % ��ֵ����ǵ������
    max_value = db(fft_max);    
    noise_value = db((sum(sum_fft(1:index - not_noise_num)) + sum(sum_fft(index + not_noise_num:end)))/(3*Nf/4));
    snr = max_value - noise_value;
end







