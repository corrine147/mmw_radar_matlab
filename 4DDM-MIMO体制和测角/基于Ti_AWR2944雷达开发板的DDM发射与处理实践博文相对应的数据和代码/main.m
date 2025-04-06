%%  主函数  控制处理过程  %%
%该函数要做的事情是：  对DDMA发射下的回波数据进行解析，看RD图，以及对目标角度进行基本的测量以验证解析的正确性%
%以3T4R为例(不包括那个俯仰向的发射天线)%


%Author: https://blog.csdn.net/xhblair?spm=1000.2115.3001.5343 %
%date: 20240118 %



clear all;  close all;  clc;

%---------雷达参数配置--------%
[RadarParament] = RadarParamentConfig;

%增加一个测角(水平向)的参数：
dr          = 0.5;    %*λ   
angleView   = linspace(-90,90,181)';
Array       = (0:1:11);    %3T4R
A           = exp(-1j*2*pi*dr*sind(angleView).*Array);   %导向矢量

%---------数据文件路径--------%
datafolder  = 'D:\Matlab仿真\4DDM-MIMO体制和测角\基于Ti_AWR2944雷达开发板的DDM发射与处理实践博文相对应的数据和代码';
binfileName_TDM = 'adc_data_TDM3T4R_Test1_Raw_0.bin';  binfilePath_TDM = strcat(datafolder,'\',binfileName_TDM);
binfileName_DDM = 'adc_data_DDM3T4R_Test1_Raw_0.bin';  binfilePath_DDM = strcat(datafolder,'\',binfileName_DDM);



%---------目标参数-----------%
% Rtarget     = 2.7;  %m  设置目标静止在雷达的正前方(当然，更标准的是我们要用CFAR找到目标，但是这里主要是确定DDMA的发射以及数据解析的正确性)
% RIndex      = floor(Rtarget/RadarParament.Rres);


%-------一帧一帧处理-------%
for frameIdx = 1:RadarParament.Numframe 

    %----------数据解析--------%
    [ADCdata_DDM] = DataParsing(binfilePath_DDM,RadarParament,frameIdx,0);   %拿出该帧的数据 大小为：距离采样点数*速度采样点数*numRx*1
    [Rangelen,dopperlen_DDM,~] = size(ADCdata_DDM);
    [ADCdata_TDM_tmp] = DataParsing(binfilePath_TDM,RadarParament,frameIdx,1);   %拿出该帧的数据 大小为：距离采样点数*速度采样点数*numRx*numTx
    %为方便后续处理，这里将之转成一个3D矩阵，其第三维的排布顺序按照虚拟阵列排布
    ADCdata_TDM(:,:,1:4) = ADCdata_TDM_tmp(:,:,1:4,1);  ADCdata_TDM(:,:,5:8) = ADCdata_TDM_tmp(:,:,1:4,2);  ADCdata_TDM(:,:,9:12) = ADCdata_TDM_tmp(:,:,1:4,3);

    %选其中一个通道的其中一个chirp的回波 看看时域数据
    figure(1);
    subplot(121);plot(1:Rangelen,real(ADCdata_DDM(:,10,2))); xlabel('采样点');ylabel('幅值');title('DDM发射模式下的时域数据');grid on;
    subplot(122);plot(1:Rangelen,real(ADCdata_TDM(:,10,2))); xlabel('采样点');ylabel('幅值');title('TDM发射模式下的时域数据');grid on;



    %----------距离FFT----------%
    RangeFFTout_DDM = fft(ADCdata_DDM,[],1);   %没有加窗
    RangeFFTout_DDM = RangeFFTout_DDM(1:round(Rangelen/2),:,:);    %因为是实采样，只要前面一半的数据
    RangeFFTout_TDM = fft(ADCdata_TDM,[],1);   %没有加窗
    RangeFFTout_TDM = RangeFFTout_TDM(1:round(Rangelen/2),:,:);    %因为是实采样，只要前面一半的数据

    %看看距离压缩后的数据
    figure(2);
    subplot(121); plot(1:round(Rangelen/2),abs(RangeFFTout_DDM(:,10,2))); xlabel('距离索引');ylabel('幅值');title('DDM发射模式下距离FFT的结果');
    subplot(122); plot(1:round(Rangelen/2),abs(RangeFFTout_TDM(:,10,2))); xlabel('距离索引');ylabel('幅值');title('TDM发射模式下距离FFT的结果');
    


    %--------速度FFT---------%
    DopplerFFTout_DDM = fftshift(fft(RangeFFTout_DDM,[],2),2);      %并做shift，把零频放在中间
%     DopplerFFTout_DDM = (fft(RangeFFTout_DDM,[],2)); 
    DopplerFFTout_TDM = fftshift(fft(RangeFFTout_TDM,[],2),2);      

    %--------DDM解调---------%   把二维FFT后的结果进行拆分：从256*256*4  变成：256*64*16   （严格来说这不是解调..）
    dopplerLength = dopperlen_DDM/4;  %4 = 发射天线个数+ emptyband数量
    DopplerFFTout_DDM_Demodulated(:,:,1:4)   = DopplerFFTout_DDM(:,1:dopplerLength,:);
    DopplerFFTout_DDM_Demodulated(:,:,13:16) = DopplerFFTout_DDM(:,dopplerLength+1:2*dopplerLength,:);
    DopplerFFTout_DDM_Demodulated(:,:,5:8)   = DopplerFFTout_DDM(:,2*dopplerLength+1:3*dopplerLength,:);
    DopplerFFTout_DDM_Demodulated(:,:,9:12)  = DopplerFFTout_DDM(:,3*dopplerLength+1:4*dopplerLength,:);


    %-------非相干积累------%
    RDmat_DDM_beforeDemodulate   = sum(abs(DopplerFFTout_DDM),3);
    RDmat_TDM   = sum(abs(DopplerFFTout_TDM),3);
    RDmat_DDM_Demodulated        = fftshift(sum(abs(DopplerFFTout_DDM_Demodulated),3),2);


    figure(3);
    subplot(231);imagesc(RDmat_TDM);xlabel('速度维');ylabel('距离维');title('TDM下二维压缩后的RD图');
    subplot(234);meshz(RDmat_TDM);xlabel('速度维');ylabel('距离维');title('TDM下二维压缩后的RD图');
    subplot(232);imagesc(RDmat_DDM_beforeDemodulate);xlabel('速度维');ylabel('距离维');title('DDM未解调前二维压缩后的RD图');
    subplot(235);meshz(RDmat_DDM_beforeDemodulate);xlabel('速度维');ylabel('距离维');title('DDM未解调前二维压缩后的RD图');
    subplot(233);imagesc(RDmat_DDM_Demodulated);xlabel('速度维');ylabel('距离维');title('DDM二维压缩后的RD图');
    subplot(236);meshz(RDmat_DDM_Demodulated);xlabel('速度维');ylabel('距离维');title('DDM二维压缩后的RD图');



    %--------测角----------%
    %拿出目标所在位置的各个点用DBF测角(这里只拿出水平维度的)
    %怎么拿出来以及排序正确很重要(同时要结合各个天线相位的设置情况)。
    %注意：我们有4个接收通道，需要从每个接收通道的矩阵中拿出NTx个点，然后进行排布，且单个接收通道中拿出来的Ntx个点其对应的位置关系是相差1个波长的。

    DopplerIndex_DDM    = [129 1 193];    %这里的顺序需要和后面的排布对应起来！
    DopplerIndex_TDM    = 33;
    RIndex              = 14;      %直接从RD图上看

    AngleData_TDM = squeeze(DopplerFFTout_TDM(RIndex,DopplerIndex_TDM,:));  %TDM12个虚拟通道的数据

    AngleData_DDM(1:4)  = DopplerFFTout_DDM(RIndex,DopplerIndex_DDM(1),:);   %因为四个接收通道间隔0.5*λ，所以第一个发射天线对应的四个通道就是找到索引后四个通道对应的数据。
    AngleData_DDM(5:8)  = DopplerFFTout_DDM(RIndex,DopplerIndex_DDM(2),:);
    AngleData_DDM(9:12) = DopplerFFTout_DDM(RIndex,DopplerIndex_DDM(3),:);
    
    P_TDM = db(abs(A*AngleData_TDM)./max(abs(A*AngleData_TDM)));
    P_DDM = db(abs(A*AngleData_DDM.')./max(abs(A*AngleData_DDM.')));


    figure(4);
    plot(angleView,P_DDM);xlabel('角度/°');ylabel('amplitude/dB');title('测角结果');grid on;hold on;
    plot(angleView,P_TDM); hold off; legend('DDM发射模式下的测角结果','TDM发射模式下的测角结果');
    
    breakpoint = 1;

end



