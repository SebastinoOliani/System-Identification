clear all
close all
clc

%% generate data
Ts = 0.1;
z = tf('z', Ts);
A = 1 + 1.1 * z^-1 + 0.7 * z^-2;
B = 0.3 * z^-1 - 1.5 * z^-2;
N = 100;

% generate dataset
rng(2);
e1 = randn(N,1) * sqrt(0.5);
rng(1);
u1 = randn(N,1) * sqrt(0.5);

y1 = lsim(B/A, u1) + e1;

figure(1)
plot(1:N, y1, 'b', 'LineWidth', 1.5)
title('Output')
xlim([1 N])
xlabel('Timestep')
ylabel('Output')

%% LS estimation of A and B
na = [1 2 3 2 3]'; % orders of A
nb = [2 2 2 3 3]'; % orders of B

Theta_LS = LS_estimation(u1, y1, na, nb);

%% instrumental variables
% generate a new dataset
rng(3);
u2 = randn(N,1) * sqrt(0.5);
rng(4);
e2 = randn(N,1) * sqrt(0.5);

% compute instrumental variable for each LS
x = zeros(N, length(na));
for i = 1:length(na)
    x(:, i) = lsim(tf([0 Theta_LS(na(i)+1:(na(i)+nb(i)), i)'], [1 Theta_LS(1:na(i), i)'], Ts), u2);
end

y2 = lsim(B/A, u2) + e2;

% estimate parameters with IV
theta_IV = zeros(na(end) + nb(end), length(na));
theta_IV = LS_estimation_IV(u2, y2, x, na, nb);

%% MSE of the output predicted by each model
% define new dataset
rng(5);
u3 = randn(N,1) * sqrt(0.5);
rng(6);
e3 = randn(N,1) * sqrt(0.5);

y3 = lsim(B/A, u3) + e3;

% predictions with LS
y_LS = zeros(N, length(na));
for i = 1:length(na)
    y_LS(:, i) = lsim(tf([0 Theta_LS(na(i)+1:(na(i)+nb(i)), i)'], [1 Theta_LS(1:na(i), i)'], Ts), u3);
end

% predictions with IV
y_IV = zeros(N, length(na));

for i = 1:length(na)
    y_IV(:, i) = lsim(tf([0 theta_IV(na(i)+1:(na(i)+nb(i)), i)'], [1 theta_IV(1:na(i), i)'], Ts), u3);
end

figure(2)
subplot(5,1,1)
plot(1:N, y3)
hold on;
plot(1:N, y_LS(:,1))
plot(1:N, y_IV(:,1))
title('n = 1, m = 2')
hold off;
legend('True plant', 'LS', 'IV')

subplot(5,1,2)
plot(1:N, y3)
hold on;
plot(1:N, y_LS(:,2))
plot(1:N, y_IV(:,2))
hold off;
title('n = 2, m = 2')

subplot(5,1,3)
plot(1:N, y3)
hold on;
plot(1:N, y_LS(:,3))
plot(1:N, y_IV(:,3))
hold off;
title('n = 3, m = 2')

subplot(5,1,4)
plot(1:N, y3)
hold on;
plot(1:N, y_LS(:,4))
plot(1:N, y_IV(:,4))
hold off;
title('n = 2, m = 3')

subplot(5,1,5)
plot(1:N, y3)
hold on;
plot(1:N, y_LS(:,5))
plot(1:N, y_IV(:,5))
hold off;
title('n = 3, m = 3')

xlabel('Timestep')

%% MSE
MSE_LS = zeros(length(na), 1);
for k = 1:length(na)
    for i = 1:N
        MSE_LS(k) = MSE_LS(k) + (y3(i) - y_LS(i, k)) ^ 2;
    end
end
MSE_LS = MSE_LS ./ N;

MSE_IV = zeros(length(na), 1);
for k = 1:length(na)
    for i = 1:N
        MSE_IV(k) = MSE_IV(k) + (y3(i) - y_IV(i, k)) ^ 2;
    end
end
MSE_IV = MSE_IV ./ N;

figure(3)
plot(1:length(na), MSE_LS, 'o')
hold on;
plot(1:length(na), MSE_IV, 'o')
hold off;
title('MSE')
xlabel('Experiment number')
legend('LS MSE', 'IV MSE')

[min_LS, idx_LS] = min(MSE_LS, [], 'all');
[min_IV, idx_IV] = min(MSE_IV, [], 'all');

if min_LS < min_IV
    disp('LS has a lower error!')
else
    disp('IV has a lower error!')
end

%% Which model has least MSE? Why?
disp('Better results with same order of the polynomials.')