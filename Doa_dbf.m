% �����������γɲ��
function angle = Doa_dbf(angle_data,array_num)
    %�������к��ź�
    array_uni = 0:1:array_num-1;                % ͬ���׾��µľ�������
    d = 0.5;                                    % ���о��ȼ��dΪ�벨��
    %ʹ��dbfɨ��
    thetascan = linspace(-90,90,1024);
    a_uni  = exp(-1i*2*pi*d*sind(thetascan)'*array_uni);           
    p_uni = angle_data.'*a_uni.'; 
    p_uni = 20*log10(abs(p_uni)./max(abs(p_uni)));
%     figure;
%     plot(thetascan,p_uni);grid on;

    [~,index] = max(p_uni);
    angle = thetascan(index);
end