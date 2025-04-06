%% ���ݽ������룬��ɴ�bin�ļ���ADC 3D���ݾ���(.mat)��ת�������������Ĳ������֤������׼ȷ�� %%
%����������������ݴ�СΪ��  numSamplePerChirp*numLoops*numRX*1  ����Ϊ����������ͬʱ���䣩

function [adcDataComplex] = DataParsing(fileFullPath, RadarParament,frameIdx,flag)

if flag == 1    %��ʾΪTDM�����ݵĽ���
   numLoops        = RadarParament.NumLoop_TDM;    
   numChirpPerLoop = RadarParament.NumTx_TDM;
else
   numLoops        = RadarParament.NumLoop_DDM;    %DDMA�¾���1���ˡ�(��Ч��һ�����߷���) 
   numChirpPerLoop = RadarParament.NumTx_DDM;  
end

numSamplePerChirp = RadarParament.NumSample;
numRX             = RadarParament.NumRx;
numIQ             = RadarParament.NumIQ;    %IQ��·�������ǵ�·��
numBits           = RadarParament.NumBits;


Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRX*numIQ;  %�����ϵ�ȫ���ĸ���

fp       = fopen(fileFullPath, 'r');
fseek(fp,(frameIdx-1)*Expected_Num_SamplesPerFrame*2, 'bof');    %������ʼ�㶨λ
DataTmp  = fread(fp,Expected_Num_SamplesPerFrame,'uint16');    %�ó���ô��������

%��ʽ����
neg                = logical(bitget(DataTmp, numBits));       %���λΪ����λ
DataTmp(neg)       = DataTmp(neg) - 2^(numBits);                     %��λΪ1��ʾΪ��������ʱ��Ҫ��ȥ2^16;

%���ݾ��������  ֻ��һ·(TDM��Ӧ����ֱ�ӿ���reshape?)
%{
%��1�� (���������⣬2944�����ݴ���������lane�ڴ���Ӧ������һ��chirp���ȴ�ǰ����Rx�������ڴ�������Rx������)
Databuffer = DataTmp;   
%Ȼ��ֱ������������
adcDataComplex     = reshape(Databuffer, numRX, numSamplePerChirp, numChirpPerLoop, numLoops);
%ֱ��reshape�����Ž���ͨ�������ŵ�chirp�µĲ������������ŷ���ͨ������������ٶ�ά��
%}

%��2����������ķ�ʽ�Ƿ���ȷ�� 
%{
Databuffer0 = reshape(DataTmp,numSamplePerChirp*4,length(DataTmp)/(numSamplePerChirp*4));
Databuffer1 = Databuffer0(1:numSamplePerChirp*2,:);
Databuffer2 = Databuffer0(numSamplePerChirp*2+1:end,:);
DataRx0 = Databuffer1(1:2:end,:);    DataRx1 = Databuffer1(2:2:end,:);
DataRx2 = Databuffer2(1:2:end,:);    DataRx3 = Databuffer2(2:2:end,:);
DataRx0 = reshape(DataRx0,numSamplePerChirp,numChirpPerLoop,numLoops);  DataRx0 = permute(DataRx0,[1 3 2]);
DataRx1 = reshape(DataRx1,numSamplePerChirp,numChirpPerLoop,numLoops);  DataRx1 = permute(DataRx1,[1 3 2]);
DataRx2 = reshape(DataRx2,numSamplePerChirp,numChirpPerLoop,numLoops);  DataRx2 = permute(DataRx2,[1 3 2]);
DataRx3 = reshape(DataRx3,numSamplePerChirp,numChirpPerLoop,numLoops);  DataRx3 = permute(DataRx3,[1 3 2]);

adcDataComplex(:,:,1,:) = DataRx0;
adcDataComplex(:,:,2,:) = DataRx1;
adcDataComplex(:,:,3,:) = DataRx2;
adcDataComplex(:,:,4,:) = DataRx3;
%}


%AWR2944Ӧ�÷��ϣ�xWR16xx/IWR6843 Real Data Format Using DCA1000
%��3�� �ȴ���һ��Rx��һ��chirp�µ��������ǵڶ���RX
datatmp2       =  reshape(DataTmp,numSamplePerChirp,numRX,numChirpPerLoop,numLoops);
adcDataComplex =  permute(datatmp2,[1 4 2 3]); 


fclose(fp);


end



























