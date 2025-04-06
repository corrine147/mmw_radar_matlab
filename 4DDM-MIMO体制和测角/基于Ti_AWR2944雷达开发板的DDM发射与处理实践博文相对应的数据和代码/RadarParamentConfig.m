%% �״�������ã����Խ���ԭʼ���� %%
%DDMA

function [RadarParament] = RadarParamentConfig


RadarParament.Numframe     = 5;     %�������֡������Ҫ���н�������֡��
RadarParament.NumTx_TDM    = 3;      %������DDMA�µȼ���ֻ��һ�����������ˣ��������ݽ���ʱ��Ҫע��
RadarParament.NumTx_DDM    = 1;  
RadarParament.NumRx        = 4;
RadarParament.NumLoop_TDM  = 64;      %��֡��ÿ���������߷����chirp��
RadarParament.NumLoop_DDM  = 256;     %��֡��ÿ���������߷����chirp��
RadarParament.NumBits      = 16;     %ADC���ݵ�λ��
RadarParament.NumIQ        = 1;      %�Ƿ�ΪIQ��·������


RadarParament.chirpSlop    = 30e12;  %chirpб��
RadarParament.fadc         = 10e6;   %Hz  ADC�Ĳ�����
RadarParament.NumSample    = 256;    %��chirp�µĲ�������
RadarParament.Buse         = RadarParament.NumSample/RadarParament.fadc*RadarParament.chirpSlop;   %ADC�������ܴ���

RadarParament.ChirpT       = 40e-6;  %chirp�ķ���ʱ��
RadarParament.Idletime     = 10e-6;  %���ʱ��
RadarParament.Tc           = RadarParament.ChirpT + RadarParament.Idletime; %ǰ���������chirp֮���ʱ�������ֵ���Է��������������ٶ�ά�ȵĲ������ڡ�
RadarParament.Frametime_DDM    = RadarParament.Tc*RadarParament.NumLoop_DDM;    %DDMA����ʱ����֡����ʱ����
RadarParament.Frametime_TDM    = RadarParament.Tc*RadarParament.NumLoop_TDM*RadarParament.NumTx_TDM; 


c                          = 3e8;
RadarParament.fc           = 77.6e9; %Hz
RadarParament.lambda       = c/RadarParament.fc;


RadarParament.Rres         = c/2/RadarParament.Buse;   %�״����۵ľ���ֱ���
FIf                        = 20e6;   %�״����Ƶ����
if RadarParament.NumIQ == 2
    Fuse = min(FIf,RadarParament.fadc);
else
    Fuse = min(FIf,RadarParament.fadc/2);
end
RadarParament.Rmax         = c*Fuse/2/RadarParament.chirpSlop;   %�����ݴ���˶��Ե�����������


RadarParament.Vres_DDM         = RadarParament.lambda/2/RadarParament.Frametime_DDM; 
RadarParament.Vres_TDM         = RadarParament.lambda/2/RadarParament.Frametime_TDM; 


%DDMA�£������ģ�����ٷ�Χ����ԭ����1/N    N�������λ�ֳɼ��ȷ֣������ʹ��empty band���ٶȽ�ģ������N = numTx
%���ʹ�ã���Ҫ���ӡ�
N = 4;   %���ֵ��Ҫ�Լ��Ӿ���������ġ�
RadarParament.Vmax_TDM         = RadarParament.lambda/4/(RadarParament.Tc*RadarParament.NumTx_TDM);
RadarParament.Vmax_DDM         = RadarParament.lambda/4/(RadarParament.Tc)/N;



end













