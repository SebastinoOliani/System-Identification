function [Theta] = LS_estimation(u, y, na, nb)

Theta = zeros(na(end) + nb(end), length(na));

for i = 1:length(na)
    row1 = zeros(na(i), 1);
    column1 = [0; -y(1:end-1)];
    row2 = zeros(nb(i), 1);
    column2 = [0; u(1:end-1)];
    
    PhiTyu = [toeplitz(column1, row1), toeplitz(column2, row2)];
    
    Theta(1:(na(i)+nb(i)), i) = PhiTyu \ y;
end

end

