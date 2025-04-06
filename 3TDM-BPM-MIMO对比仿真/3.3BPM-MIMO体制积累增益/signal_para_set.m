% �����״��źŵĲ������������ò���
function [radar_parameter] = signal_para_set()
sample_num = 384;               % ��chirp��������
slow_num = 512;                 % chirp����
c = 3e8;                        % ����
fc = 76e9;                      % ��ʼƵ��
K  = 30e12;                     % ��Ƶб��
Tp = 20e-6;                     % ��chirpʱ��
B = K*Tp;                       % �ܵ�Ƶ����
fs = 30e6;                      % ��Ƶ����Ƶ��
ts = 1/fs;                      % ����ʱ����
Ts = ts*sample_num;             % ����ʱ��,TsӦ��С��Tp
real_B = Ts*K;                  % ʵ�ʴ���
lambda = c/(fc+real_B/2);       % ����Ƶ�ʼ��㲨��
if Ts > 0.8*Tp
    printf('para has problem!');
end
% �ṹ���������
radar_parameter.c = c;
radar_parameter.Tp = Tp;
radar_parameter.fc = fc;
radar_parameter.K = K;
radar_parameter.fs = fs;
radar_parameter.ts = ts;
radar_parameter.lambda = lambda;
radar_parameter.sample_num = sample_num;
radar_parameter.slow_num = slow_num;
radar_parameter.real_B = real_B;
end