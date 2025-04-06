%% 写一个函数实现：
%输入：Tx个数、emptyband的个数、emptyband的位置、需要发射的chirp个数。
%输出：每个Tx在360°范围内应该设置的理想相位值、以及可以设置的实际相位值(因为相位值的步进值是5.626°)

clear; close all; clc;

%参数设置
PhaseMin      = 5.625;
numTx         = 3;
numChirps     = 256;

numEmptybands = 1;

EmptybandIndex = 1;   
%为1表示处在Tx1与Tx2之间，为2表示处在Tx2与Tx3之间，以此类推。 (不过无论几个emptyband，这里默认它们是连续的)
%不过不管emptyband是如何放置的，后续具体使用时找到Tx对应的行即可。



%相位值计算
N_all     =  numTx + numEmptybands;
PhaseSetp =  360/N_all;

Idea_Phase = zeros(N_all,numChirps);   %预设一个矩阵装载相位。
Real_Phase = zeros(N_all,numChirps);   %实际的相位   
IndexUse   = zeros(N_all,numChirps);   %实际相位对应的应该设置的index值： 从0--63

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



%画图看看结果  但是跳过emptyband
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
hold off;grid on; legend('Tx1','Tx2','Tx3'); %需要更改

figure(2); hold on;
ii = 1;
while ii <= N_all
    if ii == IndexSkip(1)
        ii = IndexSkip(end)+1;
    end
    plot(1:numChirps,IndexUse(ii,:));
    ii = ii+1;
end  
hold off;grid on;legend('Tx1','Tx2','Tx3'); %需要更改




















