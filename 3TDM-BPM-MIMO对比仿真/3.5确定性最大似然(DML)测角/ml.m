% 最大似然估计算法仿真似然函数和对数似然函数的谱图
% Demonstration of Maximum Likelihood Estimation in Matlab
%   Author: Mathuranathan (https://www.gaussianwaves.com)
%   License : creative commons : Attribution-NonCommercial-ShareAlike 3.0 Unported
clc;
close all;
clear;
 
A=1.3;
N=50; %采样点
x=A+randn(1,N);
s=1; %标准差 s=1 
 
rangeA=-2:0.0001:5; %搜索范围
L=zeros(1,length(rangeA)); %搜索长度
 
for i=1:length(rangeA)
    %计算似然函数
    L(i) = exp(-sum((x-rangeA(i)).^2)/(2*s^2));  %忽略常数项
end
 
[maxL,index]=max(L); %选择极大似然函数的最大值
display('Maximum Likelihood of A');
display(rangeA(index));
 
%Plotting Commands
figure;
plot(rangeA,L);hold on;
stem(rangeA(index),L(index),'r'); %Point the Maximum Likelihood Estimate
displayText=['\leftarrow Likelihood of A=' num2str(rangeA(index))];
title('Maximum Likelihood Estimation of Parameter A');
xlabel('\leftarrow A');
ylabel('Likelihood');
text(rangeA(index),L(index)/3,displayText,'HorizontalAlignment','left');
 
figure;
plot(rangeA,log(L));hold on;
YL = ylim;YMIN = YL(1);
plot([rangeA(index) rangeA(index)],[YMIN log(L(index))] ,'r'); %Point the Maximum Likelihood Estimate
title('Log Likelihood Function');
xlabel('\leftarrow A');
ylabel('Log Likelihood');
text([rangeA(index)],YMIN/2,displayText,'HorizontalAlignment','left');