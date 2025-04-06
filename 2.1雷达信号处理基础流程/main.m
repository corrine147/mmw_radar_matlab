% FMCWĿ��ز���FFT��FFT_1D��FFT_2D��CFAR����
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
% tar_inf = [15 -2;
%            105 3;
%            70.8 20;
%            18 -1;
%            90.3 -15;];          % ����Ŀ�����,�ٶ�
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

ADC_data = awgn(ADC_data,10);        % ��������,0dB
% ADC_data = real(ADC_data);        % ֻʹ���ź�ʵ��
% figure;
% plot(real(ADC_data(:,1)));
% figure;
% plot(imag(ADC_data(:,1)));

% 1D FFT
FFT_1D = fft(ADC_data,sample_num,1);
R_label = R_res*(1:sample_num);
[XX,YY] = meshgrid(slow_time*1e3,R_label);
figure;plot(R_label,db(abs(FFT_1D(:,1))));hold on;
plot([tar_inf(:,1),tar_inf(:,1)],ylim,'m-.');legend('����ֲ�','��ʵ����');
title('����ֲ�');xlabel('Range/m');ylabel('power(dB)');grid on;
figure;mesh(XX,YY,db(abs(FFT_1D)));colorbar;title('һάFFT');xlabel('time/ms');ylabel('range/m');zlabel('power(dB)');grid on;
figure;imagesc(slow_time*1e3,R_label,(abs(FFT_1D)));colorbar;title('һάFFT');xlabel('time/ms');ylabel('range/m');zlabel('power(dB)');

% 2D FFT
FFT_2D = fftshift(fft(FFT_1D,[],2),2);
V_lable = (-V_Max:V_res:(V_Max-V_res));
[XX,YY] = meshgrid(V_lable,R_label);
figure;mesh(XX,YY,db(abs(FFT_2D)));colorbar;title('��άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;imagesc(V_lable,R_label,(abs(FFT_2D)));colorbar;title('��άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');

% CFAR �����ټ�
result_thre1 = zeros(sample_num,slow_num);
echo_thre = zeros(sample_num,slow_num);
for j = 1:slow_num
    echo_thre(:,j) = CFAR(FFT_2D(:,j),16);
    result = abs(FFT_2D(:,j)) - echo_thre(:,j);
    [loc,~] = find(result>1);
    result_thre1(loc,j) = 1;
end
result_thre2 = zeros(sample_num,slow_num);
for j = 1:sample_num
    echo_thre(j,:) = CFAR(FFT_2D(j,:),12);
    result = abs(FFT_2D(j,:)) - echo_thre(j,:);
    [~,loc] = find(result>1);
    result_thre2(j,loc) = 1;
end
result_thre3 = result_thre1 .* result_thre2;
figure;mesh(XX,YY,result_thre3);colorbar;title('��άFFT + CFAR');xlabel('speed(m/s)');ylabel('range/m');zlabel('amplitude');
figure;imagesc(V_lable,R_label,result_thre3);colorbar;title('��άFFT + CFAR');xlabel('speed(m/s)');ylabel('range/m');zlabel('amplitude');
% CFAR������Ա�
[m,n] = find(result_thre3>0.8);
target_v = V_lable(n);
target_r = R_label(m);
figure;plot(tar_inf(:,2),tar_inf(:,1),'r*');hold on;
plot(target_v,target_r,'d');legend('��ʵĿ������ٶ�','���Ŀ������ٶ�');grid on;title('CFAR������Ա�');



