% 多通道幅相一致性仿真，相位不一致和校准
clc;clear;close all;
N = 16;                 % 通道数
L = 1000;               % 时域信号序列长度      
fs = 60000;             % 采样频率
ts = 1/fs;              % 采样时间间隔
t = ts:ts:L*ts;         % 时间序列
phase = 0:0.1*pi:(N-1)*0.1*pi;  % 每个通道的初始相位，单位弧度
A = 2.5:2.5:N*2.5;
s = zeros(N, L);        
ss = zeros(N, L);       % 增加直流分量用于区分显示
f = 259;                % 信号频率
 
for i = 1:N
    s(i,:) = cos(2*pi*f*t + phase(i)) + 0.1*randn(1,L);
    ss(i,:) = A(i) + s(i,:);
end
% 画出时域信号
figure;
plot(t, s);
grid on;
title('混叠在一起的16通道信号');
figure;
plot(t, ss);
grid on;
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16');
title('用直流幅度区分的16通道信号');
 
Nf = 65536;              % 计算FFT的点数
df = fs/Nf;              % 频率分辨率
s_fft = zeros(N,Nf);     % s的FFT结果 
for i = 1:N
    s_fft(i,:) = 10*log10(abs(fft(s(i,:),Nf)));
end
% 画出FFT幅值结果，这里频率太低FFT结果不明显
% figure;
% for i = 1:N
%     subplot(4,4,i);
%     plot((1:Nf/8)*df,s_fft(i,(1:Nf/8)));
%     ylim([-10,30]);
%     grid on;
%     title(num2str(i));
% end
% suptitle('FFT频率分析-dB');
 
% 通过最大幅值响应计算频率
[fft_max, index] = max(s_fft,[],2);
f_cal = index*df;
f_cal = mean(f_cal);
f_error = f_cal - f;
 
% 相移精度，单位°
d_phase = 360*f*ts;
min_dots = 360 / d_phase;
corr_num = floor(min_dots) + 1;
s_corr_sum = zeros(N,corr_num);
% 通过自相关（互相关）计算相移，这里只取一个周期用来计算
for i = 1:N
    for j = 1:corr_num
        s_corr = s(1,j:j+corr_num-1).*s(i,1:corr_num);
        s_corr_sum(i,j) = sum(s_corr);
    end
end
% 画出时域信号相关曲线
figure;
for i = 1:N
    subplot(4,4,i);
    plot(1:corr_num,s_corr_sum(i,:));
    grid on;
    title(num2str(i));
end
suptitle('第1通道为基准的相关曲线');
 
% 取自相关最大值索引为相移值
[max_sum,index] = max(s_corr_sum,[],2);
phase_cal = index*d_phase;
phase_error = phase_cal - phase'*180/pi;      % 这里单位转成°
phase_mean_error = mean(phase_error(2:N));

% 补偿后的信号
s_comp = zeros(N, L);  
ss_comp = zeros(N, L);
for i = 1:N
    s_comp(i,:) = cos(2*pi*f*t + phase_error(i)*pi/180) + 0.1*randn(1,L);
    ss_comp(i,:) = A(i) + s_comp(i,:);
end
% 画出补偿后的时域信号
figure;
plot(t,s_comp);
grid on;
title('混叠在一起的补偿后16通道信号');
figure;
plot(t,ss_comp);grid on;
title('各通道补偿后的时域曲线叠加');
grid on;
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16');
title('用直流幅度区分的补偿后16通道信号');




