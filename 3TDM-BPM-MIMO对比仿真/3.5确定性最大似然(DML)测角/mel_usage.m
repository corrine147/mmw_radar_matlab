
clear;close all;clc;
x = randn(10000,1);
h = histogram(x);   % 分布直方图；统计某个数在整个数的占比是多少，区间为坐标左闭右开的统计
phat = mle(x);


