% 均匀线阵波束形成测角
function angle = Doa_dbf(angle_data,array_num)
    %构造阵列和信号
    array_uni = 0:1:array_num-1;                % 同样孔径下的均匀阵列
    d = 0.5;                                    % 阵列均匀间隔d为半波长
    %使用dbf扫描
    thetascan = linspace(-90,90,1024);
    a_uni  = exp(-1i*2*pi*d*sind(thetascan)'*array_uni);           
    p_uni = angle_data.'*a_uni.'; 
    p_uni = 20*log10(abs(p_uni)./max(abs(p_uni)));
%     figure;
%     plot(thetascan,p_uni);grid on;

    [~,index] = max(p_uni);
    angle = thetascan(index);
end