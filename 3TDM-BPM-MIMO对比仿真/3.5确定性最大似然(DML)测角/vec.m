% 将矩阵转换成向量，按列排序
function ret = vec(X)
    [m,n] = size(X);
    ret = zeros(1,m*n);
    for i = 1:m
        for j = 1:n
            ret((i-1)*n + j) = X(i,j);
        end
    end
end
