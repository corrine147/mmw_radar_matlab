
clear;close all;clc;
x = randn(10000,1);
h = histogram(x);   % �ֲ�ֱ��ͼ��ͳ��ĳ��������������ռ���Ƕ��٣�����Ϊ��������ҿ���ͳ��
phat = mle(x);


