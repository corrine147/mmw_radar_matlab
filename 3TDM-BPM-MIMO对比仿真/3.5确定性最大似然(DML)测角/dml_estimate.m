% ȷ���������Ȼ���Ʒ���
clc; clear; close all;

%% ��������
%%% ����Ƶ��
lambda = 1;         % ����
k = 2*pi/lambda;    % ����
%%% ���в���
N = 10;                 % ��Ԫ����
d = 0.5*lambda;         % ��Ԫ��� 
z = (0:d:(N-1)*d)';     % ��Ԫ����ֲ�
n_list = (0:N-1)';
%%% �ź�Դ����
% target_theta = [-30, -10.6, 20.4, 60.2]';     % �����������
target_theta = [-30, 35]';     % �����������
M = length(target_theta);           % �ź�Դ��Ŀ
%%% �������
SNR = 10;     % �����(dB)
L = 500;     % ��������

%% ���н����źŷ���ģ��
A = exp(1j*k*z*sind(target_theta'));          % ���;���
S = randn(M, L);         % �����ź�
X = A*S;                        % ���н����ź�
X1 = awgn(X, SNR, 'measured');      % ���ظ�˹������������н����ź�
% ���н����źŵ�Э�������
R = X1*X1'/L;    

%% ȷ����ʼֵ
theta_list = (-90:0.1:90)'*pi/180;
for idx = 1:length(theta_list)
    %%% ����DML�׺���
    psi0 = k*d*sin(theta_list(idx));        % electrical angle
    A = exp(1j*n_list*psi0');
    P_A = A*inv(A'*A)*A';                       % projection matrix wrt V
    f(idx) = trace(P_A*R);
end
%%% �׺�����һ��
f = real(f);
f = (f-min(f))/(max(f)-min(f));
[~, ii] = max(f);       % ��ֵ
theta0 = theta_list(ii);
psi = k*d*sin(theta0);      % ����ֵ��Ϊţ�ٵ�����ʼֵ
%%% ��ȡ�׷�ֵ
[f_peaks, f_peaks_idx] = findpeaks(f);     % ��ȡ��ֵ
[f_peaks, I] = sort(f_peaks, 'descend');    % ��ֵ��������
f_peaks_idx = f_peaks_idx(I);
f_peaks = f_peaks(1:M);             % ��ȡǰM��
f_peaks_idx = f_peaks_idx(1:M);
theta_est = theta_list(f_peaks_idx)*180/pi;   % ���Ʒ���
disp('�ź�Դ���Ʒ���Ϊ��');
disp(theta_est);
%%% ���� DML �׺���
figure;
plot(theta_list*180/pi, db(f));hold on;
xlabel('\theta (deg)');ylabel('DML�ռ���');title('ȷ���������Ȼ���ƽǶ�');
plot([target_theta,target_theta]',ylim,'m-.');hold on;
plot(theta_est, db(f_peaks), 'r.', 'MarkerSize', 15);hold on;
for idx = 1:M
    text(theta_est(idx)+3, f_peaks(idx), sprintf('%0.1f��', theta_est(idx)));
end
grid on;ylim([-90,10]);