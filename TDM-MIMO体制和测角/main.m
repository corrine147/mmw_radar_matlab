% FMCW目标回波、FFT、FFT_1D、FFT_2D、CFAR仿真
clc;clear;close all;
% 目标参数
tar_inf = [15 1 10;
           18 -5 0;
           70.8 20 0;
           90.3 -15 20;
           105 10 -30;];          % 仿真目标距离,速度,角度
% 雷达参数设置
radar_parameter = signal_para_set();
c = radar_parameter.c;
k = radar_parameter.K;
Tp = radar_parameter.Tp;
fs = radar_parameter.fs;
lambda = radar_parameter.lambda;
sample_num = radar_parameter.sample_num;        	% 单chirp采样点数
slow_num = radar_parameter.slow_num;                % chirp数量

R_res = c/2/radar_parameter.real_B;             % 距离分辨率
R_Max = R_res*sample_num;                       % 最远距离
V_Max = radar_parameter.lambda/4/(2*Tp);        % 最大不模糊速度,相比于单发模式降低了一半
V_res = V_Max*2/slow_num;                       % 速度分辨率，相比于单发模式也降低了一半
% 这里仅考虑均匀排布天线阵列
Tx_num = 2;                                     % 发射天线数
Rx_num = 4;                                     % 接收天线数
dRx = lambda/2;
Channel_data = zeros(Tx_num*Rx_num,sample_num,slow_num);
for Tx_id = 1:Tx_num
    for Rx_id = 1:Rx_num
        ADC_data = zeros(sample_num,slow_num);
        for temp_tar = 1:size(tar_inf,1)
            slow_time = (Tx_id - 1 + 2*(0:slow_num-1))*Tp;
            tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;        % 生成波形，距离R=R0+v*dt
            % temp_data = signal_generate(1, 2*tar_pos, radar_parameter);
            fc_array = radar_parameter.fc*ones(sample_num,1);
            fast_time = ((0:sample_num-1)/fs)';
            delay = 2*tar_pos/c;
            phi = ((Tx_id-1)*Rx_num+Rx_id)*dRx/lambda * sind(tar_inf(temp_tar,3));
            temp_data = exp(1i*2*pi*(phi + k*fast_time*delay + fc_array*delay  - 1/2*k*(delay).^2));
            ADC_data = ADC_data + temp_data;
        end
        ADC_data = awgn(ADC_data,0);        % 产生噪声,0dB
        % 1D FFT
        FFT_1D = fft(ADC_data,sample_num,1);
        % 2D FFT
        FFT_2D = fftshift(fft(FFT_1D,[],2),2);
        Channel_data((Tx_id-1)*Rx_num+Rx_id,:,:) = FFT_2D;
    end
end

% 显示各个通道的RD图
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
Acc_channel_data = zeros(sample_num,slow_num);
figure;
for i = 1:Tx_num*Rx_num
    subplot(2,4,i);
    temp_data = abs(squeeze(Channel_data(i,:,:)));
    Acc_channel_data = Acc_channel_data + temp_data;
    imagesc(V_label,R_label,temp_data);title('单通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
end

[XX,YY] = meshgrid(V_label,R_label);
figure;mesh(XX,YY,db(temp_data));colorbar;title('第8通道二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;mesh(XX,YY,db(Acc_channel_data));colorbar;title('积累后二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;imagesc(V_label,R_label,Acc_channel_data);colorbar;title('积累后二维FFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');

% CFAR检测点迹
result_thre1 = zeros(sample_num,slow_num);
echo_thre = zeros(sample_num,slow_num);
for j = 1:slow_num
    echo_thre(:,j) = CFAR(Acc_channel_data(:,j),16);
    result = abs(Acc_channel_data(:,j)) - echo_thre(:,j);
    [loc,~] = find(result>1);
    result_thre1(loc,j) = 1;
end
result_thre2 = zeros(sample_num,slow_num);
for j = 1:sample_num
    echo_thre(j,:) = CFAR(Acc_channel_data(j,:),12);
    result = abs(Acc_channel_data(j,:)) - echo_thre(j,:);
    [~,loc] = find(result>1);
    result_thre2(j,loc) = 1;
end
result_thre3 = result_thre1 .* result_thre2;
figure;mesh(XX,YY,result_thre3);colorbar;title('积累后二维FFT + CFAR');xlabel('speed(m/s)');ylabel('range/m');zlabel('amplitude');
figure;imagesc(V_label,R_label,result_thre3);colorbar;title('积累后二维FFT + CFAR');xlabel('speed(m/s)');ylabel('range/m');zlabel('amplitude');
% CFAR检测结果RD对比
[m,n] = find(result_thre3>0.8);
target_r = R_label(m);
target_v = V_label(n);
figure;plot(tar_inf(:,2),tar_inf(:,1),'r*');hold on;xlim([-V_Max,V_Max]);ylim([0,R_Max]);
plot(target_v,target_r,'d');legend('真实目标距离速度','检测目标距离速度');grid on;title('CFAR检测结果RD对比');
% CFAR检测结果RA比
target_num = length(target_r);
angle_data = zeros(1, Tx_num*Rx_num);
target_a = zeros(1,target_num);                     % 未补偿的目标角度
Comp_phi = zeros(1, Tx_num*Rx_num);                 % 补偿的相位
target_a_comp = zeros(1,target_num);                % 补偿后的目标角度
for i = 1:target_num
    angle_data = Channel_data(:,m(i),n(i));
    target_a(i) = Doa_dbf(angle_data,Tx_num*Rx_num);
    % 计算速度引起的相位偏差补偿
    phi = 4*pi*V_label(n(i))*radar_parameter.fc*Tp/c;
    for Tx_id = 1:Tx_num
        for Rx_id = 1:Rx_num
            Comp_phi((Tx_id-1)*Rx_num + Rx_id) = exp(-1i*(Tx_id - 1)*phi);
        end
    end
    angle_data = angle_data.*Comp_phi.';
    % 计算补偿后的角度
    target_a_comp(i) = Doa_dbf(angle_data,Tx_num*Rx_num);
end
figure;plot(tar_inf(:,3),tar_inf(:,1),'r*');hold on;xlim([-60,60]);ylim([0,R_Max]);
plot(target_a,target_r,'d');legend('真实目标距离角度','检测目标距离角度');grid on;title('CFAR检测结果RA对比');
figure;plot(tar_inf(:,3),tar_inf(:,1),'r*');hold on;xlim([-60,60]);ylim([0,R_Max]);
plot(target_a_comp,target_r,'go');legend('真实目标距离角度','检测目标距离角度');grid on;title('相位补偿后CFAR检测结果RA对比');








