function [u] = generate_u(e_u, N)
% initialization
u_1 = 0; u_2 = 0;
e_1 = 0; e_2 = 0;
u = zeros(N,1);

for k = 1:N
    u(k) = 0.1 * u_1 + 0.12 * u_2 + e_1 + 0.2 * e_2;
    u_2 = u_1; u_1 = u(k);
    e_2 = e_1; e_1 = e_u(k);
end

end

