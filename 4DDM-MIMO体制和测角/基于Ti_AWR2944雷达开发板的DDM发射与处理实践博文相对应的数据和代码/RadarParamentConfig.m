%% 雷达参数设置，用以解析原始数据 %%
%DDMA

function [RadarParament] = RadarParamentConfig


RadarParament.Numframe     = 5;     %发射的总帧数或想要进行解析的总帧数
RadarParament.NumTx_TDM    = 3;      %不过在DDMA下等价于只有一个发射天线了，这在数据解析时需要注意
RadarParament.NumTx_DDM    = 1;  
RadarParament.NumRx        = 4;
RadarParament.NumLoop_TDM  = 64;      %单帧下每个发射天线发射的chirp数
RadarParament.NumLoop_DDM  = 256;     %单帧下每个发射天线发射的chirp数
RadarParament.NumBits      = 16;     %ADC数据的位数
RadarParament.NumIQ        = 1;      %是否为IQ两路采样？


RadarParament.chirpSlop    = 30e12;  %chirp斜率
RadarParament.fadc         = 10e6;   %Hz  ADC的采样率
RadarParament.NumSample    = 256;    %单chirp下的采样点数
RadarParament.Buse         = RadarParament.NumSample/RadarParament.fadc*RadarParament.chirpSlop;   %ADC采样的总带宽

RadarParament.ChirpT       = 40e-6;  %chirp的发射时长
RadarParament.Idletime     = 10e-6;  %间隔时间
RadarParament.Tc           = RadarParament.ChirpT + RadarParament.Idletime; %前后发射的两个chirp之间的时长，这个值乘以发射天线数才是速度维度的采样周期。
RadarParament.Frametime_DDM    = RadarParament.Tc*RadarParament.NumLoop_DDM;    %DDMA发射时，单帧发射时长。
RadarParament.Frametime_TDM    = RadarParament.Tc*RadarParament.NumLoop_TDM*RadarParament.NumTx_TDM; 


c                          = 3e8;
RadarParament.fc           = 77.6e9; %Hz
RadarParament.lambda       = c/RadarParament.fc;


RadarParament.Rres         = c/2/RadarParament.Buse;   %雷达理论的距离分辨率
FIf                        = 20e6;   %雷达的中频带宽
if RadarParament.NumIQ == 2
    Fuse = min(FIf,RadarParament.fadc);
else
    Fuse = min(FIf,RadarParament.fadc/2);
end
RadarParament.Rmax         = c*Fuse/2/RadarParament.chirpSlop;   %从数据处理端而言的最大测量距离


RadarParament.Vres_DDM         = RadarParament.lambda/2/RadarParament.Frametime_DDM; 
RadarParament.Vres_TDM         = RadarParament.lambda/2/RadarParament.Frametime_TDM; 


%DDMA下，最大无模糊测速范围会变成原来的1/N    N是你把相位分成几等分？如果不使用empty band做速度解模糊，则N = numTx
%如果使用，则要增加。
N = 4;   %这个值需要自己视具体情况更改。
RadarParament.Vmax_TDM         = RadarParament.lambda/4/(RadarParament.Tc*RadarParament.NumTx_TDM);
RadarParament.Vmax_DDM         = RadarParament.lambda/4/(RadarParament.Tc)/N;



end













