%%  ������  ���ƴ������  %%
%�ú���Ҫ���������ǣ�  ��DDMA�����µĻز����ݽ��н�������RDͼ���Լ���Ŀ��ǶȽ��л����Ĳ�������֤��������ȷ��%
%��3T4RΪ��(�������Ǹ�������ķ�������)%


%Author: https://blog.csdn.net/xhblair?spm=1000.2115.3001.5343 %
%date: 20240118 %



clear all;  close all;  clc;

%---------�״��������--------%
[RadarParament] = RadarParamentConfig;

%����һ�����(ˮƽ��)�Ĳ�����
dr          = 0.5;    %*��   
angleView   = linspace(-90,90,181)';
Array       = (0:1:11);    %3T4R
A           = exp(-1j*2*pi*dr*sind(angleView).*Array);   %����ʸ��

%---------�����ļ�·��--------%
datafolder  = 'D:\Matlab����\4DDM-MIMO���ƺͲ��\����Ti_AWR2944�״￪�����DDM�����봦��ʵ���������Ӧ�����ݺʹ���';
binfileName_TDM = 'adc_data_TDM3T4R_Test1_Raw_0.bin';  binfilePath_TDM = strcat(datafolder,'\',binfileName_TDM);
binfileName_DDM = 'adc_data_DDM3T4R_Test1_Raw_0.bin';  binfilePath_DDM = strcat(datafolder,'\',binfileName_DDM);



%---------Ŀ�����-----------%
% Rtarget     = 2.7;  %m  ����Ŀ�꾲ֹ���״����ǰ��(��Ȼ������׼��������Ҫ��CFAR�ҵ�Ŀ�꣬����������Ҫ��ȷ��DDMA�ķ����Լ����ݽ�������ȷ��)
% RIndex      = floor(Rtarget/RadarParament.Rres);


%-------һ֡һ֡����-------%
for frameIdx = 1:RadarParament.Numframe 

    %----------���ݽ���--------%
    [ADCdata_DDM] = DataParsing(binfilePath_DDM,RadarParament,frameIdx,0);   %�ó���֡������ ��СΪ�������������*�ٶȲ�������*numRx*1
    [Rangelen,dopperlen_DDM,~] = size(ADCdata_DDM);
    [ADCdata_TDM_tmp] = DataParsing(binfilePath_TDM,RadarParament,frameIdx,1);   %�ó���֡������ ��СΪ�������������*�ٶȲ�������*numRx*numTx
    %Ϊ��������������ｫ֮ת��һ��3D���������ά���Ų�˳�������������Ų�
    ADCdata_TDM(:,:,1:4) = ADCdata_TDM_tmp(:,:,1:4,1);  ADCdata_TDM(:,:,5:8) = ADCdata_TDM_tmp(:,:,1:4,2);  ADCdata_TDM(:,:,9:12) = ADCdata_TDM_tmp(:,:,1:4,3);

    %ѡ����һ��ͨ��������һ��chirp�Ļز� ����ʱ������
    figure(1);
    subplot(121);plot(1:Rangelen,real(ADCdata_DDM(:,10,2))); xlabel('������');ylabel('��ֵ');title('DDM����ģʽ�µ�ʱ������');grid on;
    subplot(122);plot(1:Rangelen,real(ADCdata_TDM(:,10,2))); xlabel('������');ylabel('��ֵ');title('TDM����ģʽ�µ�ʱ������');grid on;



    %----------����FFT----------%
    RangeFFTout_DDM = fft(ADCdata_DDM,[],1);   %û�мӴ�
    RangeFFTout_DDM = RangeFFTout_DDM(1:round(Rangelen/2),:,:);    %��Ϊ��ʵ������ֻҪǰ��һ�������
    RangeFFTout_TDM = fft(ADCdata_TDM,[],1);   %û�мӴ�
    RangeFFTout_TDM = RangeFFTout_TDM(1:round(Rangelen/2),:,:);    %��Ϊ��ʵ������ֻҪǰ��һ�������

    %��������ѹ���������
    figure(2);
    subplot(121); plot(1:round(Rangelen/2),abs(RangeFFTout_DDM(:,10,2))); xlabel('��������');ylabel('��ֵ');title('DDM����ģʽ�¾���FFT�Ľ��');
    subplot(122); plot(1:round(Rangelen/2),abs(RangeFFTout_TDM(:,10,2))); xlabel('��������');ylabel('��ֵ');title('TDM����ģʽ�¾���FFT�Ľ��');
    


    %--------�ٶ�FFT---------%
    DopplerFFTout_DDM = fftshift(fft(RangeFFTout_DDM,[],2),2);      %����shift������Ƶ�����м�
%     DopplerFFTout_DDM = (fft(RangeFFTout_DDM,[],2)); 
    DopplerFFTout_TDM = fftshift(fft(RangeFFTout_TDM,[],2),2);      

    %--------DDM���---------%   �Ѷ�άFFT��Ľ�����в�֣���256*256*4  ��ɣ�256*64*16   ���ϸ���˵�ⲻ�ǽ��..��
    dopplerLength = dopperlen_DDM/4;  %4 = �������߸���+ emptyband����
    DopplerFFTout_DDM_Demodulated(:,:,1:4)   = DopplerFFTout_DDM(:,1:dopplerLength,:);
    DopplerFFTout_DDM_Demodulated(:,:,13:16) = DopplerFFTout_DDM(:,dopplerLength+1:2*dopplerLength,:);
    DopplerFFTout_DDM_Demodulated(:,:,5:8)   = DopplerFFTout_DDM(:,2*dopplerLength+1:3*dopplerLength,:);
    DopplerFFTout_DDM_Demodulated(:,:,9:12)  = DopplerFFTout_DDM(:,3*dopplerLength+1:4*dopplerLength,:);


    %-------����ɻ���------%
    RDmat_DDM_beforeDemodulate   = sum(abs(DopplerFFTout_DDM),3);
    RDmat_TDM   = sum(abs(DopplerFFTout_TDM),3);
    RDmat_DDM_Demodulated        = fftshift(sum(abs(DopplerFFTout_DDM_Demodulated),3),2);


    figure(3);
    subplot(231);imagesc(RDmat_TDM);xlabel('�ٶ�ά');ylabel('����ά');title('TDM�¶�άѹ�����RDͼ');
    subplot(234);meshz(RDmat_TDM);xlabel('�ٶ�ά');ylabel('����ά');title('TDM�¶�άѹ�����RDͼ');
    subplot(232);imagesc(RDmat_DDM_beforeDemodulate);xlabel('�ٶ�ά');ylabel('����ά');title('DDMδ���ǰ��άѹ�����RDͼ');
    subplot(235);meshz(RDmat_DDM_beforeDemodulate);xlabel('�ٶ�ά');ylabel('����ά');title('DDMδ���ǰ��άѹ�����RDͼ');
    subplot(233);imagesc(RDmat_DDM_Demodulated);xlabel('�ٶ�ά');ylabel('����ά');title('DDM��άѹ�����RDͼ');
    subplot(236);meshz(RDmat_DDM_Demodulated);xlabel('�ٶ�ά');ylabel('����ά');title('DDM��άѹ�����RDͼ');



    %--------���----------%
    %�ó�Ŀ������λ�õĸ�������DBF���(����ֻ�ó�ˮƽά�ȵ�)
    %��ô�ó����Լ�������ȷ����Ҫ(ͬʱҪ��ϸ���������λ���������)��
    %ע�⣺������4������ͨ������Ҫ��ÿ������ͨ���ľ������ó�NTx���㣬Ȼ������Ų����ҵ�������ͨ�����ó�����Ntx�������Ӧ��λ�ù�ϵ�����1�������ġ�

    DopplerIndex_DDM    = [129 1 193];    %�����˳����Ҫ�ͺ�����Ų���Ӧ������
    DopplerIndex_TDM    = 33;
    RIndex              = 14;      %ֱ�Ӵ�RDͼ�Ͽ�

    AngleData_TDM = squeeze(DopplerFFTout_TDM(RIndex,DopplerIndex_TDM,:));  %TDM12������ͨ��������

    AngleData_DDM(1:4)  = DopplerFFTout_DDM(RIndex,DopplerIndex_DDM(1),:);   %��Ϊ�ĸ�����ͨ�����0.5*�ˣ����Ե�һ���������߶�Ӧ���ĸ�ͨ�������ҵ��������ĸ�ͨ����Ӧ�����ݡ�
    AngleData_DDM(5:8)  = DopplerFFTout_DDM(RIndex,DopplerIndex_DDM(2),:);
    AngleData_DDM(9:12) = DopplerFFTout_DDM(RIndex,DopplerIndex_DDM(3),:);
    
    P_TDM = db(abs(A*AngleData_TDM)./max(abs(A*AngleData_TDM)));
    P_DDM = db(abs(A*AngleData_DDM.')./max(abs(A*AngleData_DDM.')));


    figure(4);
    plot(angleView,P_DDM);xlabel('�Ƕ�/��');ylabel('amplitude/dB');title('��ǽ��');grid on;hold on;
    plot(angleView,P_TDM); hold off; legend('DDM����ģʽ�µĲ�ǽ��','TDM����ģʽ�µĲ�ǽ��');
    
    breakpoint = 1;

end



