function [Theta] = LS_estimation_IV(u, y, x, na, nb)
Theta = zeros(na(end) + nb(end), length(na));

for i = 1:length(na)
    row1 = zeros(na(i),1);
    row2 = zeros(nb(i), 1);
    column = [0; -x(1:end-1, i)];
    input = [0; u(1:end-1)];
    Z = [toeplitz(column, row1), toeplitz(input, row2)];

    column1 = [0; -y(1:end-1)];
    
    PhiTyu = [toeplitz(column1, row1), toeplitz(input, row2)];
    
    Theta(1:na(i)+nb(i), i) = (Z' * PhiTyu) \ Z' * y;
end

end

