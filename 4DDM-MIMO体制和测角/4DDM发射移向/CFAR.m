%% CA-CFAR
function [Result_Map] = CFAR(Echo_PC,Ref_Num)
Pro_Num = 7;
% Ref_Num = 2^4;                                    % 可修改
P_fa = 5e-11;                                       % 虚警率
Threshold = Ref_Num*(P_fa^(-1/(2*Ref_Num))-1);      % 理论门限计算
Num_All = length(Echo_PC);
for i_Num = 1:Num_All
    if i_Num<=Pro_Num+1
        Result_Map(i_Num) = mean(abs(Echo_PC(i_Num+Pro_Num+1:i_Num+Pro_Num+Ref_Num)));
    elseif Pro_Num+1<i_Num && i_Num<=(Pro_Num+Ref_Num)
        Result_Map_Left = mean(abs(Echo_PC(1:i_Num-Pro_Num-1)));
        Result_Map_Right = mean(abs(Echo_PC(i_Num+Pro_Num+1:i_Num+Pro_Num+Ref_Num)));
        Result_Map(i_Num) = mean([Result_Map_Left,Result_Map_Right]);
    elseif (Pro_Num+Ref_Num)<i_Num && i_Num<=(Num_All-Pro_Num-Ref_Num)
        Result_Map_Left = mean(abs(Echo_PC(i_Num-Pro_Num-Ref_Num:i_Num-Pro_Num-1)));
        Result_Map_Right = mean(abs(Echo_PC(i_Num+Pro_Num+1:i_Num+Pro_Num+Ref_Num)));
        Result_Map(i_Num) = mean([Result_Map_Left,Result_Map_Right]);
    elseif (Num_All-Pro_Num-Ref_Num)<i_Num && i_Num<(Num_All-Pro_Num)
        Result_Map_Left = mean(abs(Echo_PC(i_Num-Pro_Num-Ref_Num:i_Num-Pro_Num-1)));
        Result_Map_Right = mean(abs(Echo_PC(i_Num+Pro_Num+1:end)));
        Result_Map(i_Num) = mean([Result_Map_Left,Result_Map_Right]);
    elseif i_Num>=(Num_All-Pro_Num)
        Result_Map(i_Num) = mean(abs(Echo_PC(i_Num-Pro_Num-Ref_Num:i_Num-Pro_Num-1)));
    end
    X_Ref = [];
    Result_Map_Left = [];
    Result_Map_Right = [];
end
Result_Map = Result_Map.*Threshold;
end