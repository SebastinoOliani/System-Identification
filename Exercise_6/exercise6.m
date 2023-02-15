clear all
clc

%%
a=5/8;
b=11/10;
cov=0.05;
t_max=80;
% number of period
n=2;
% PRBS in range -2;2
u = generate_u(n ,t_max);

%% 1: generate output signal
v1 = randn(n*t_max,1)*sqrt(cov);
y1 = generate_y(u, v1, a, b, n, t_max);
g1 = zeros(t_max,1);
c = u(1:n*t_max,1);
r = [u(1), zeros(1, t_max-1)];
phi = toeplitz(c,r);
g1 = (phi'*phi)\phi'*y1;

v2 = randn(n*t_max,1)*sqrt(cov);
y2 = generate_y(u, v2, a, b, n, t_max);
g2 = zeros(t_max,1);
c = u(1:n*t_max,1);
r = [u(1), zeros(1, t_max-1)];
phi = toeplitz(c,r);
g2 = (phi'*phi)\phi'*y2;

%% computing g as average
g = 0.5*(g1+g2);

%% 2: apply N periods of the signal
g = zeros(t_max,29);
g_0 = zeros(t_max,29);
for N=2:30
    v = randn(N*t_max,1)*sqrt(cov);
    u = generate_u(N, t_max);
    y = generate_y(u, v, a, b, N, t_max);
    y_0 = generate_y(u, zeros(N*t_max,1), a, b, N, t_max);
    c = u(1:N*t_max,1);
    r = [u(1), zeros(1, t_max-1)];
    phi = toeplitz(c,r);
    g(:,N-1) = (phi'*phi)\phi'*y;
    g_0(:, N-1) = (phi'*phi)\phi'*y_0;
end

%% plot error gN - g0 as function of N
figure(1)
error = zeros(29,1);
average = ones(29,1);

for N=2:30
    error(N-1) = norm(g(:,N-1)-g_0(:,N-1));
    sum = 0;
    for i = 1:t_max
        sum = sum + error(N-1);
    end
    if N==30
        average = average * sum/N;
    end
end

N=2:30;
plot(N, error)
hold on;
plot(N, average)
% title('Error')
% xlabel('N, number of periods')
% ylabel('Error')
% xlim([N(1), N(29)])
% legend('Period', 'Average')
% hold off;

%% 3: input generated as rand
g = zeros(t_max,29);
g_0 = zeros(t_max,29);
for N=2:30
    v = randn(N*t_max,1)*sqrt(cov);
    u = rand(N*t_max);
    y = generate_y(u, v, a, b, N, t_max);
    y_0 = generate_y(u, zeros(N*t_max,1), a, b, N, t_max);
    c = u(1:N*t_max,1);
    r = [u(1), zeros(1, t_max-1)];
    phi = toeplitz(c,r);
    g(:,N-1) = (phi'*phi)\phi'*y;
    g_0(:, N-1) = (phi'*phi)\phi'*y_0;
end

% figure(2)
error = zeros(29,1);
average = ones(29,1);

for N=2:30
    error(N-1) = norm(g(:,N-1)-g_0(:,N-1));
    sum = 0;
    for i = 1:t_max
        sum = sum + error(N-1);
    end
    if N==30
        average = average * sum/N;
    end
end

N=2:30;
plot(N, error)
% hold on;
plot(N, average)
title('Error')
xlabel('N, number of periods')
ylabel('Error')
xlim([N(1), N(29)])
legend('Period', 'Average', 'Period rand', 'Average rand')
hold off;
% uniform random noise increases the error

%% 4: repeat M times part 1
M = 100;
g = zeros(t_max,M);
g_M = zeros(t_max,M);
for i=1:M
    v = randn(n*t_max,1)*sqrt(cov);
    y = generate_y(u, v, a, b, n, t_max);
    c = u(1:n*t_max,1);
    r = [u(1), zeros(1, t_max-1)];
    phi = toeplitz(c,r);
    g(:,i) = (phi'*phi)\phi'*y;
    for j=1:i
        g_M(:,i) = g_M(:,i) + g_M(:,j);
    end
    g_M(:,i) = g_M(:,i) + g(:,i);
end
g_M = g_M ./ M;

% generate g_0
y_0 = generate_y(u, zeros(n*t_max,1), a, b, n, t_max);
g_0 = (phi'*phi)\phi'*y_0;

error = zeros(M,1);
for i=1:M
    error(i) = norm(g_M(:,i)-g_0(:));
%     sum = 0;
%     for i = 1:t_max
%         sum = sum + error(N-1);
%     end
%     if N==30
%         average = average * sum/N;
%     end
end

figure(2)
% plot with axis in logaritmic scale to see the rate of growth
semilogy(1:M, error)
xlabel('M')
ylabel('Error')
title('Error')
xlim([1,M])
% the bigger the period, the bigger the error