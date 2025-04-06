% ��ͨ������һ���Է��棬��λ��һ�º�У׼
clc;clear;close all;
N = 16;                 % ͨ����
L = 1000;               % ʱ���ź����г���      
fs = 60000;             % ����Ƶ��
ts = 1/fs;              % ����ʱ����
t = ts:ts:L*ts;         % ʱ������
phase = 0:0.1*pi:(N-1)*0.1*pi;  % ÿ��ͨ���ĳ�ʼ��λ����λ����
A = 2.5:2.5:N*2.5;
s = zeros(N, L);        
ss = zeros(N, L);       % ����ֱ����������������ʾ
f = 259;                % �ź�Ƶ��
 
for i = 1:N
    s(i,:) = cos(2*pi*f*t + phase(i)) + 0.1*randn(1,L);
    ss(i,:) = A(i) + s(i,:);
end
% ����ʱ���ź�
figure;
plot(t, s);
grid on;
title('�����һ���16ͨ���ź�');
figure;
plot(t, ss);
grid on;
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16');
title('��ֱ���������ֵ�16ͨ���ź�');
 
Nf = 65536;              % ����FFT�ĵ���
df = fs/Nf;              % Ƶ�ʷֱ���
s_fft = zeros(N,Nf);     % s��FFT��� 
for i = 1:N
    s_fft(i,:) = 10*log10(abs(fft(s(i,:),Nf)));
end
% ����FFT��ֵ���������Ƶ��̫��FFT���������
% figure;
% for i = 1:N
%     subplot(4,4,i);
%     plot((1:Nf/8)*df,s_fft(i,(1:Nf/8)));
%     ylim([-10,30]);
%     grid on;
%     title(num2str(i));
% end
% suptitle('FFTƵ�ʷ���-dB');
 
% ͨ������ֵ��Ӧ����Ƶ��
[fft_max, index] = max(s_fft,[],2);
f_cal = index*df;
f_cal = mean(f_cal);
f_error = f_cal - f;
 
% ���ƾ��ȣ���λ��
d_phase = 360*f*ts;
min_dots = 360 / d_phase;
corr_num = floor(min_dots) + 1;
s_corr_sum = zeros(N,corr_num);
% ͨ������أ�����أ��������ƣ�����ֻȡһ��������������
for i = 1:N
    for j = 1:corr_num
        s_corr = s(1,j:j+corr_num-1).*s(i,1:corr_num);
        s_corr_sum(i,j) = sum(s_corr);
    end
end
% ����ʱ���ź��������
figure;
for i = 1:N
    subplot(4,4,i);
    plot(1:corr_num,s_corr_sum(i,:));
    grid on;
    title(num2str(i));
end
suptitle('��1ͨ��Ϊ��׼���������');
 
% ȡ��������ֵ����Ϊ����ֵ
[max_sum,index] = max(s_corr_sum,[],2);
phase_cal = index*d_phase;
phase_error = phase_cal - phase'*180/pi;      % ���ﵥλת�ɡ�
phase_mean_error = mean(phase_error(2:N));

% ��������ź�
s_comp = zeros(N, L);  
ss_comp = zeros(N, L);
for i = 1:N
    s_comp(i,:) = cos(2*pi*f*t + phase_error(i)*pi/180) + 0.1*randn(1,L);
    ss_comp(i,:) = A(i) + s_comp(i,:);
end
% �����������ʱ���ź�
figure;
plot(t,s_comp);
grid on;
title('�����һ��Ĳ�����16ͨ���ź�');
figure;
plot(t,ss_comp);grid on;
title('��ͨ���������ʱ�����ߵ���');
grid on;
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16');
title('��ֱ���������ֵĲ�����16ͨ���ź�');




