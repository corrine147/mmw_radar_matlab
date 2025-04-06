% FMCWĿ��ز���FFT��FFT_1D��FFT_2D��CFAR����
clc;clear;close all;
% Ŀ�����
tar_inf = [15 1 10;
           18 -5 0;
           70.8 20 0;
           90.3 -15 20;
           105 10 -30;];          % ����Ŀ�����,�ٶ�,�Ƕ�
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
V_Max = radar_parameter.lambda/4/(2*Tp);        % ���ģ���ٶ�,����ڵ���ģʽ������һ��
V_res = V_Max*2/slow_num;                       % �ٶȷֱ��ʣ�����ڵ���ģʽҲ������һ��
% ��������Ǿ����Ų���������
Tx_num = 2;                                     % ����������
Rx_num = 4;                                     % ����������
dRx = lambda/2;
Channel_data = zeros(Tx_num*Rx_num,sample_num,slow_num);
for Tx_id = 1:Tx_num
    for Rx_id = 1:Rx_num
        ADC_data = zeros(sample_num,slow_num);
        for temp_tar = 1:size(tar_inf,1)
            slow_time = (Tx_id - 1 + 2*(0:slow_num-1))*Tp;
            tar_pos = tar_inf(temp_tar,1) + tar_inf(temp_tar,2)*slow_time;        % ���ɲ��Σ�����R=R0+v*dt
            % temp_data = signal_generate(1, 2*tar_pos, radar_parameter);
            fc_array = radar_parameter.fc*ones(sample_num,1);
            fast_time = ((0:sample_num-1)/fs)';
            delay = 2*tar_pos/c;
            phi = ((Tx_id-1)*Rx_num+Rx_id)*dRx/lambda * sind(tar_inf(temp_tar,3));
            temp_data = exp(1i*2*pi*(phi + k*fast_time*delay + fc_array*delay  - 1/2*k*(delay).^2));
            ADC_data = ADC_data + temp_data;
        end
        ADC_data = awgn(ADC_data,0);        % ��������,0dB
        % 1D FFT
        FFT_1D = fft(ADC_data,sample_num,1);
        % 2D FFT
        FFT_2D = fftshift(fft(FFT_1D,[],2),2);
        Channel_data((Tx_id-1)*Rx_num+Rx_id,:,:) = FFT_2D;
    end
end

% ��ʾ����ͨ����RDͼ
V_label = (-V_Max:V_res:(V_Max-V_res));
R_label = R_res*(1:sample_num);
Acc_channel_data = zeros(sample_num,slow_num);
figure;
for i = 1:Tx_num*Rx_num
    subplot(2,4,i);
    temp_data = abs(squeeze(Channel_data(i,:,:)));
    Acc_channel_data = Acc_channel_data + temp_data;
    imagesc(V_label,R_label,temp_data);title('��ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');
end

[XX,YY] = meshgrid(V_label,R_label);
figure;mesh(XX,YY,db(temp_data));colorbar;title('��8ͨ����άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;mesh(XX,YY,db(Acc_channel_data));colorbar;title('���ۺ��άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');grid on;
figure;imagesc(V_label,R_label,Acc_channel_data);colorbar;title('���ۺ��άFFT');xlabel('speed(m/s)');ylabel('range/m');zlabel('power(dB)');

% CFAR���㼣
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
figure;mesh(XX,YY,result_thre3);colorbar;title('���ۺ��άFFT + CFAR');xlabel('speed(m/s)');ylabel('range/m');zlabel('amplitude');
figure;imagesc(V_label,R_label,result_thre3);colorbar;title('���ۺ��άFFT + CFAR');xlabel('speed(m/s)');ylabel('range/m');zlabel('amplitude');
% CFAR�����RD�Ա�
[m,n] = find(result_thre3>0.8);
target_r = R_label(m);
target_v = V_label(n);
figure;plot(tar_inf(:,2),tar_inf(:,1),'r*');hold on;xlim([-V_Max,V_Max]);ylim([0,R_Max]);
plot(target_v,target_r,'d');legend('��ʵĿ������ٶ�','���Ŀ������ٶ�');grid on;title('CFAR�����RD�Ա�');
% CFAR�����RA��
target_num = length(target_r);
angle_data = zeros(1, Tx_num*Rx_num);
target_a = zeros(1,target_num);                     % δ������Ŀ��Ƕ�
Comp_phi = zeros(1, Tx_num*Rx_num);                 % ��������λ
target_a_comp = zeros(1,target_num);                % �������Ŀ��Ƕ�
for i = 1:target_num
    angle_data = Channel_data(:,m(i),n(i));
    target_a(i) = Doa_dbf(angle_data,Tx_num*Rx_num);
    % �����ٶ��������λƫ���
    phi = 4*pi*V_label(n(i))*radar_parameter.fc*Tp/c;
    for Tx_id = 1:Tx_num
        for Rx_id = 1:Rx_num
            Comp_phi((Tx_id-1)*Rx_num + Rx_id) = exp(-1i*(Tx_id - 1)*phi);
        end
    end
    angle_data = angle_data.*Comp_phi.';
    % ���㲹����ĽǶ�
    target_a_comp(i) = Doa_dbf(angle_data,Tx_num*Rx_num);
end
figure;plot(tar_inf(:,3),tar_inf(:,1),'r*');hold on;xlim([-60,60]);ylim([0,R_Max]);
plot(target_a,target_r,'d');legend('��ʵĿ�����Ƕ�','���Ŀ�����Ƕ�');grid on;title('CFAR�����RA�Ա�');
figure;plot(tar_inf(:,3),tar_inf(:,1),'r*');hold on;xlim([-60,60]);ylim([0,R_Max]);
plot(target_a_comp,target_r,'go');legend('��ʵĿ�����Ƕ�','���Ŀ�����Ƕ�');grid on;title('��λ������CFAR�����RA�Ա�');








