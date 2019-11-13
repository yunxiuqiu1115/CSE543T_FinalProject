% kernel matrix

function k = kernelMatrix(a, b, kernel)
    n = length(a);
    m = length(b);
    k = zeros(n,m);
    for j = 1:m
        for i = 1:n
            k(i,j) = kernel(a(i), b(j));
        end
    end
end