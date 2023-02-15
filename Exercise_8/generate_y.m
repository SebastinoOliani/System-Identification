function [y] = generate_y(u, e, N, A, B, C)
% initialization
y = zeros(N,1);
u_1 = 0; u_2 = 0;
e_1 = 0; e_2 = 0;
y_1 = 0; y_2 = 0;

for k = 1:N
    y(k) = -A(2) * y_1 - A(3) * y_2 + B(2) * u_1 + B(3) * u_2 + ...
        e(k) + C(2) * e_1 + C(3) * e_2;
    u_2 = u_1; u_1 = u(k);
    e_2 = e_1; e_1 = e(k);
    y_2 = y_1; y_1 = y(k);
end

end

