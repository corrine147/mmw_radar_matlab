% ʹ��ţ�ٷ��������Ȼ����
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
target_theta = [-30, -6.6, 20.4, 60.2]';     % �����������
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
theta_list = (-90:1:90)'*pi/180;
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
%%% ��ȡ�׷�ֵ
[f_peaks, f_peaks_idx] = findpeaks(f);     % ��ȡ��ֵ
[f_peaks, I] = sort(f_peaks, 'descend');    % ��ֵ��������
f_peaks_idx = f_peaks_idx(I);
f_peaks = f_peaks(1:M);             % ��ȡǰ M ��
f_peaks_idx = f_peaks_idx(1:M);
theta0 = theta_list(f_peaks_idx);  
psi = k*d*sin(theta0);      % ����ֵ��Ϊţ�ٵ�����ʼֵ

%% ȷ���������Ȼ�㷨DML + ţ�ٷ� ��Ŀ�ģ�ʵ�ָ��߷ֱ棩
%%% ������������
delta = 1e-4;       % ������ֹ����
del_psi = 1e3;
%%% ţ�ٷ�����
idx = 1;
while del_psi > delta
    %%% �������о���
    A = exp(1j*n_list*psi');
    %%% �������о����α�������ͶӰ����
    A_plus = inv(A'*A)*A';
    P_A1 = eye(N) - A*A_plus;
    %%% �������о��� A ��һ�׵����� D
    D = diag(1j*n_list)*A;
    %%% һ�׵���
    F = -2*real(diag(A_plus*R*P_A1*D));
    %%% Hessian ����
    H = 2*real((D'*P_A1*D).*((A_plus*R*A_plus').'));
    %%% ���´��Ż����� \psi   
    del_psi = inv(H)*F;
    psi = psi - del_psi;
    %%% ����
    del_psi = abs(del_psi);
    idx = idx + 1;
end
theta_est = asin(psi/(k*d))*180/pi;     % �Ƕȹ���ֵ
fprintf('ţ�ٷ������ź�Դ����Ϊ��');
disp(sort(theta_est)');
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