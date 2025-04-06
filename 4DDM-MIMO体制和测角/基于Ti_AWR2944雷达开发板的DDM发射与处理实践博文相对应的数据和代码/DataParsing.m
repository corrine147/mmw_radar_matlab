%% 数据解析代码，完成从bin文件到ADC 3D数据矩阵(.mat)的转换，并做基本的测距以验证解析的准确性 %%
%这里解析出来的数据大小为：  numSamplePerChirp*numLoops*numRX*1  （因为是所有天线同时发射）

function [adcDataComplex] = DataParsing(fileFullPath, RadarParament,frameIdx,flag)

if flag == 1    %表示为TDM下数据的解析
   numLoops        = RadarParament.NumLoop_TDM;    
   numChirpPerLoop = RadarParament.NumTx_TDM;
else
   numLoops        = RadarParament.NumLoop_DDM;    %DDMA下就是1个了。(等效于一个天线发射) 
   numChirpPerLoop = RadarParament.NumTx_DDM;  
end

numSamplePerChirp = RadarParament.NumSample;
numRX             = RadarParament.NumRx;
numIQ             = RadarParament.NumIQ;    %IQ两路采样还是单路？
numBits           = RadarParament.NumBits;


Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRX*numIQ;  %理论上的全部的个数

fp       = fopen(fileFullPath, 'r');
fseek(fp,(frameIdx-1)*Expected_Num_SamplesPerFrame*2, 'bof');    %数据起始点定位
DataTmp  = fread(fp,Expected_Num_SamplesPerFrame,'uint16');    %拿出这么多数据来

%格式解析
neg                = logical(bitget(DataTmp, numBits));       %最高位为符号位
DataTmp(neg)       = DataTmp(neg) - 2^(numBits);                     %首位为1表示为负数，此时需要减去2^16;

%数据矩阵的制作  只有一路(TDM下应该是直接可以reshape?)
%{
%法1： (但是有问题，2944的数据传输是两个lane在传，应该是在一个chirp下先传前两个Rx的数据在传后两个Rx的数据)
Databuffer = DataTmp;   
%然后直接下面这样：
adcDataComplex     = reshape(Databuffer, numRX, numSamplePerChirp, numChirpPerLoop, numLoops);
%直接reshape，先排接收通道，再排单chirp下的采样点数，再排发射通道，最后再排速度维度
%}

%法2：这个解析的方式是否正确？ 
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


%AWR2944应该符合：xWR16xx/IWR6843 Real Data Format Using DCA1000
%法3： 先传完一个Rx的一个chirp下的数据再是第二个RX
datatmp2       =  reshape(DataTmp,numSamplePerChirp,numRX,numChirpPerLoop,numLoops);
adcDataComplex =  permute(datatmp2,[1 4 2 3]); 


fclose(fp);


end



























