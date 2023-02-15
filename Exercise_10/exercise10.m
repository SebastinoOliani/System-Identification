clear all
close all
clc

%% model
s = tf('s');
Gs1 = (s^2 + 2 * 0.1 * 3 * s + 3 ^ 2) / (s^2 + 2 * 0.1 * 3.5 * s + 3.5 ^ 2);
Gs2 = 5000 / ((s + 50) * (s + 200));
Gs = Gs1 * Gs2;
Ts = 0.02;

Gdz = c2d(Gs, Ts, 'zoh');

z = tf('z', Ts);
Cdz = (1.25 * z - 0.75) / (z - 1);

%% sensitivity and complementary functions
Sdz = 1 / (1 + Cdz * Gdz);
Tdz = 1 - Sdz;

%% identification with PRBS 
N = 1023; % period of the signal
n = 500; % number of periods

y_average = zeros(N, 1);
u_average = zeros(N, 1);

r_prbs = 0.1 * idinput(N, 'PRBS');
r = zeros(n * N, 1);
for k = 1:n
    r((k-1) * N + 1 :(k * N)) = r_prbs;
end

v = randn(N * n, 1) * sqrt(0.1);

y = lsim(Tdz, r) + lsim(Sdz, v);
u = lsim(Sdz * Cdz, r) - lsim(Sdz * Cdz, v);

% remove first period
new_y = y(N+1 : end);
new_u = u(N+1 : end);

for i = 0:n-2
    y_average = y_average + new_y(i * N + 1: (i+1) * N);
    u_average = u_average + new_u(i * N + 1: (i+1) * N);
end

y_average = y_average ./ (n - 1);
u_average = u_average ./ (n - 1);

%% estimate the ETFE
fs = 1:1:N*1/2;

Y = fft(y_average);
U = fft(u_average);
U = U(fs); Y = Y(fs);
omega = 2*pi/N/Ts * (1:(N-1)/2);

G = Y ./ U;

Gdz_freq_resp = freqresp(Gdz, omega);
Gdz_freq_resp = reshape(Gdz_freq_resp, length(fs), 1);

figure (1)
subplot(2,1,1)
loglog(omega, abs(Gdz_freq_resp))
hold on;
loglog(omega, abs(G))
hold off;
legend('True plant', 'ETFE estimate')
grid on;
xlim([omega(1) omega(end)])

subplot(2,1,2)
semilogx(omega, angle(Gdz_freq_resp) * 180 / pi)
hold on;
semilogx(omega, angle(G) * 180 / pi)
hold off;
grid on;
xlim([omega(1) omega(end)])

%% estimate Sdz
e_average = zeros(N, 1);

r_prbs = 0.1 * idinput(N, 'PRBS');
r = zeros(n * N, 1);
for k = 1:n
    r((k-1) * N + 1 :(k * N)) = r_prbs;
end

% noise
v = randn(N * n, 1) * sqrt(0.1);

y = lsim(Tdz, r) + lsim(Sdz, v);
e = r - y;

% remove first period
new_y = y(N+1 : end);
new_e = e(N+1 : end);

for i = 0:n-2
    y_average = y_average + new_y(i * N + 1: (i+1) * N);
    e_average = e_average + new_e(i * N + 1: (i+1) * N);
end

y_average = y_average ./ (n - 1);
e_average = e_average ./ (n - 1);

E = fft(e_average);
R = fft(r_prbs);
E = E(fs); R = R(fs);
S = E ./ R;

Sdz_freq_resp = freqresp(Sdz, omega);
Sdz_freq_resp = reshape(Sdz_freq_resp, length(fs), 1);

figure (2)
subplot(2,1,1)
loglog(omega, abs(Sdz_freq_resp))
hold on;
plot(omega, abs(S))
hold off;
legend('True sensitivity function', 'ETFE estimate')
grid on;
xlim([omega(1) omega(end)])

subplot(2,1,2)
semilogx(omega, angle(Sdz_freq_resp) * 180 / pi)
hold on;
plot(omega, angle(S) * 180 / pi)
hold off;
grid on;
xlim([omega(1) omega(end)])

%% excitation signal @ controller output
y_average = zeros(N, 1);
w_average = zeros(N, 1);

w_prbs = 0.1 * idinput(N, 'PRBS');
w = zeros(n * N, 1);
for k = 1:n
    w((k-1) * N + 1 :(k * N)) = w_prbs;
end
% noise
v = randn(N * n, 1) * sqrt(0.1);

y = lsim(Sdz * Gdz, w) + lsim(Sdz, v);

% remove first period
new_y = y(N+1 : end);
new_w = w(N+1 : end);

for i = 0:n-2
    y_average = y_average + new_y(i * N + 1: (i+1) * N);
    w_average = w_average + new_w(i * N + 1: (i+1) * N);
end

y_average = y_average ./ (n - 1);
w_average = w_average ./ (n - 1);

Y = fft(y_average);
W = fft(w_average);
Y = Y(fs); W = W(fs);

% complementary sensitivity function
Gdz_Sdz = Y ./ W;

Gdz_Sdz_freq_resp = freqresp(Sdz * Gdz, omega);
Gdz_Sdz_freq_resp = reshape(Gdz_Sdz_freq_resp, length(fs), 1);

figure (3)
subplot(2,1,1)
loglog(omega, abs(Gdz_Sdz_freq_resp))
hold on;
plot(omega, abs(Gdz_Sdz))
hold off;
legend('True complementary function', 'ETFE estimate')
grid on;
xlim([omega(1) omega(end)])

subplot(2,1,2)
semilogx(omega, angle(Gdz_Sdz_freq_resp) * 180 / pi)
hold on;
plot(omega, angle(Gdz_Sdz) * 180 / pi)
hold off;
grid on;
xlim([omega(1) omega(end)])

%% estimate Gdz
C = U ./ E;

figure (4)
subplot(2,1,1)
loglog(omega, abs(Gdz_freq_resp))
hold on;
plot(omega, abs(Gdz_Sdz ./ S))
plot(omega, abs(Gdz_Sdz ./ (1 - Gdz_Sdz .* C)))
hold off;
legend('True plant', 'ETFE estimate', 'Indirect method')
grid on;
xlim([omega(1) omega(end)])

subplot(2,1,2)
semilogx(omega, angle(Gdz_freq_resp) * 180 / pi)
hold on;
plot(omega, angle(Gdz_Sdz ./ S) * 180 / pi)
plot(omega, angle(Gdz_Sdz ./ (1 - Gdz_Sdz .* C)) * 180 / pi)
hold off;
grid on;
xlim([omega(1) omega(end)])