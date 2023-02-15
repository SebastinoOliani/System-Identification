clear all
close all
clc

%% Initialization
N = 10^4;
z = tf('z');
A_ = 1 - 1.5 * z^-1 + 0.7 * z^-2;
A = [1 -1.5 0.7]';
B_ = 1.0 * z^-1 + 0.5 * z^-2;
B = [0 1.0 0.5]';
C_ = 1 - 1.0 * z^-1 + 0.2 * z^-2;
C = [1 -1.0 0.2]';

e = randn(N, 1) * sqrt(1);
e_u = randn(N, 1) * sqrt(1);

%% Fit input sequence u
u = generate_u(e_u, N);

% compute y
% y = lsim(B_/A_, u) + lsim(C_/A_, e);
y = generate_y(u, e, N, A, B, C);

%% Asymptotically unbiased LS
% u filtered
u_f = lsim(C_^-1, u);
% y filtered
y_f = lsim(C_^-1, y);

% regressor
PhiTyu = zeros(N,4);
PhiTyu(1,:) = [0, 0, 0, 0];
PhiTyu(2,:) = [-y_f(1), 0, u_f(1), 0];
for k=3:N
    PhiTyu(k, :) = [-y_f(k-1), -y_f(k-2), u_f(k-1), u_f(k-2)];
end

Theta = PhiTyu \ y_f; % Theta = [a1 a2 b1 b2];

%% Predicted values of last 1000 points
% plot can be made with filtered or non filtered signals
LS = (Theta(3) * z^-1 + Theta(4) * z^-2)/(1 + Theta(1) * z^-1 + Theta(2) * z^-2);
y_p = lsim(LS, u_f) + lsim(C_/(1 + Theta(1) * z^-1 + Theta(2) * z^-2), e);

k = N-999:N;
plot(k, y_f(k))
hold on;
plot(k, y_p(k))
title('Predicted output')
xlabel('Timestep')
ylabel('Output')
legend('True output', 'Predicted output')
hold off;

%% 100 different realizations of e(k)
reali = 100;
theta_r = zeros(4,reali);

for i=1:reali
    e = randn(N, 1) * sqrt(1);
    e_u = randn(N, 1) * sqrt(1);
    u = generate_u(e_u, N);
    y = generate_y(u, e, N, A, B, C);
    u_f = lsim(C_^-1, u);
    y_f = lsim(C_^-1, y);

    PhiTyu = zeros(N,4);
    PhiTyu(1,:) = [0, 0, 0, 0];
    PhiTyu(2,:) = [-y_f(1), 0, u_f(1), 0];
    for k=3:N
        PhiTyu(k, :) = [-y_f(k-1), -y_f(k-2), u_f(k-1), u_f(k-2)];
    end

    theta_r(:, i) = PhiTyu \ y_f;
end

figure (2)
histogram(theta_r(1,:))
title('a1')

figure (3)
histogram(theta_r(2,:))
title('a2')

figure (4)
histogram(theta_r(3,:))
title('b1')

figure (5)
histogram(theta_r(4,:))
title('b2')