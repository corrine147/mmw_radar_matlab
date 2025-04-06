%% дһ������ʵ�֣�
%���룺Tx������emptyband�ĸ�����emptyband��λ�á���Ҫ�����chirp������
%�����ÿ��Tx��360�㷶Χ��Ӧ�����õ�������λֵ���Լ��������õ�ʵ����λֵ(��Ϊ��λֵ�Ĳ���ֵ��5.626��)

clear; close all; clc;

%��������
PhaseMin      = 5.625;
numTx         = 3;
numChirps     = 256;

numEmptybands = 1;

EmptybandIndex = 1;   
%Ϊ1��ʾ����Tx1��Tx2֮�䣬Ϊ2��ʾ����Tx2��Tx3֮�䣬�Դ����ơ� (�������ۼ���emptyband������Ĭ��������������)
%��������emptyband����η��õģ���������ʹ��ʱ�ҵ�Tx��Ӧ���м��ɡ�



%��λֵ����
N_all     =  numTx + numEmptybands;
PhaseSetp =  360/N_all;

Idea_Phase = zeros(N_all,numChirps);   %Ԥ��һ������װ����λ��
Real_Phase = zeros(N_all,numChirps);   %ʵ�ʵ���λ   
IndexUse   = zeros(N_all,numChirps);   %ʵ����λ��Ӧ��Ӧ�����õ�indexֵ�� ��0--63

for ii = 1:N_all
    PhaseSetp_use = (ii-1)*PhaseSetp;
    for jj = 1:numChirps
       Phasetmp   = rem( (0 + (jj-1)*PhaseSetp_use),360);
       Idea_Phase(ii,jj) = Phasetmp;
       Quatient  = round(Phasetmp/PhaseMin);
       Real_Phase(ii,jj) = Quatient*PhaseMin;
       IndexUse(ii,jj)   = Quatient;
    end
end



%��ͼ�������  ��������emptyband
IndexSkip = (EmptybandIndex+1 : EmptybandIndex+numEmptybands) ;
figure(1); hold on;
ii = 1;
while ii <= N_all
    if ii == IndexSkip(1)
        ii = IndexSkip(end)+1;
    end
    plot(1:numChirps,Idea_Phase(ii,:));
    plot(1:numChirps,Real_Phase(ii,:));
%     scatter(1:numChirps,Idea_Phase(ii,:));
%     scatter(1:numChirps,Real_Phase(ii,:));
    ii = ii+1;
end  
hold off;grid on; legend('Tx1','Tx2','Tx3'); %��Ҫ����

figure(2); hold on;
ii = 1;
while ii <= N_all
    if ii == IndexSkip(1)
        ii = IndexSkip(end)+1;
    end
    plot(1:numChirps,IndexUse(ii,:));
    ii = ii+1;
end  
hold off;grid on;legend('Tx1','Tx2','Tx3'); %��Ҫ����




















