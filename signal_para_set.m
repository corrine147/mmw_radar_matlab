% 设置雷达信号的参数，返回配置参数
function [radar_parameter] = signal_para_set()
sample_num = 384;               % 单chirp采样点数
slow_num = 512;                 % chirp数量
c = 3e8;                        % 光速
fc = 76e9;                      % 起始频率
K  = 30e12;                     % 调频斜率
Tp = 20e-6;                     % 单chirp时长
B = K*Tp;                       % 总调频带宽
fs = 30e6;                      % 中频采样频率
ts = 1/fs;                      % 采样时间间隔
Ts = ts*sample_num;             % 采样时间,Ts应当小于Tp
real_B = Ts*K;                  % 实际带宽
lambda = c/(fc+real_B/2);       % 中心频率计算波长
if Ts > 0.8*Tp
    printf('para has problem!');
end
% 结构体参数设置
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