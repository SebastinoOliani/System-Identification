clear all
close all
clc

%% identify the system for different ARX orders
z = tf('z');
A = 1 + 0.25 * z^-1 - 0.2 * z^-2 + 0.1 * z^-3 + 0.05 * z^-4;
B = 0.6 * z^-1 + 0.3 * z^-2 - 0.05 * z^-3;

n = 100;
sigma_u = 1;
u = randn(n,1) * sigma_u;
sigma_d = 0.1;
d = randn(n,1) * sigma_d;
d_filter = lsim(A^-1, d);
y_filter = lsim(B/A, u) + d_filter;


% initial conditions
y_id = zeros(n,1);
y_id1 = 0; y_id2 = 0; y_id3 = 0; y_id4 = 0;
u_1 = 0; u_2 = 0; u_3 = 0;

% identified input
for i=1:n
    y_id(i) = - 0.25 * y_id1 + 0.2 * y_id2 - 0.1 * y_id3 - 0.05 * y_id4...
        + 0.6 * u_1 + 0.3 * u_2 - 0.05 * u_3 + d(i);
    % update previous values
    y_id4 = y_id3; y_id3 = y_id2; y_id2 = y_id1; y_id1 = y_id(i);
    u_3 = u_2; u_2 = u_1; u_1 = u(i);
end

%%
z_id = [y_id,u]; 

for m=1:8
    sys_id{m} = arx(z_id, [m m 1]);
end

%% obtain prediction for each identified ARX models
n = 100;
sigma_u_val = 1;
u_val = randn(n, 1) * sigma_u_val;
sigma_d_val = 0.1;
d_val = randn(n, 1) * sigma_d_val;

% initial conditions 
y_val = zeros(n,1);
y_val1 = 0; y_val2 = 0; y_val3 = 0; y_val4 = 0;
u_1 = 0; u_2 = 0; u_3 = 0;

% validation output
for i=1:n
    y_val(i) = - 0.25 * y_val1 + 0.2 * y_val2 - 0.1 * y_val3 - 0.05 * y_val4...
        + 0.6 * u_1 + 0.3 * u_2 - 0.05 * u_3 + d_val(i);
    % update previous values
    y_val4 = y_val3; y_val3 = y_val2; y_val2 = y_val1; y_val1 = y_val(i);
    u_3 = u_2; u_2 = u_1; u_1 = u_val(i);
end

%%
z_val = [y_val, u_val]; 

for m=1:8
    sys_val{m} = arx(z_val, [m m 1]);
end

%%
figure (1)
Y = zeros(n,1);
A_ = [1 z^-1 z^-2 z^-3 z^-4 z^-5 z^-6 z^-7 z^-8];
B_ = [0 z^-1 z^-2 z^-3 z^-4 z^-5 z^-6 z^-7 z^-8];

% predicted outputs with identification data set
for i=1:8
    arx_model = (B_(1, 1:i+1) * sys_id{i}.B') / (A_(1, 1:i+1) * sys_id{i}.A');
    Y = lsim(arx_model, u);
    plot(1:n, Y)
    hold on;
end

% true noisy output
plot(1:n, y_id, 'o')
title('Predicted identification output')
xlabel('Timestep')
ylabel('y')
hold off;
legend('ARX 1', 'ARX 2', 'ARX 3', 'ARX 4', 'ARX 5', 'ARX 6', 'ARX 7', 'ARX 8', 'True noisy output')

%%
figure (2)
Y_val = zeros(n,1);

% predicted outputs with validation data set
for i=1:8
   arx_model = (B_(1, 1:i+1) * sys_val{i}.B') / (A_(1, 1:i+1) * sys_val{i}.A');
   Y_val = lsim(arx_model, u_val);
   plot(1:n, Y_val)
   hold on;
end

% true noisy output
plot(1:n, y_val, 'o')
title('Predicted validation output')
xlabel('Timestep')
ylabel('y')
hold off;
legend('ARX 1', 'ARX 2', 'ARX 3', 'ARX 4', 'ARX 5', 'ARX 6', 'ARX 7', 'ARX 8', 'True noisy output')

%% compute mean square prediction error
figure (3)
% identification data set
MSE_id=zeros(8,1);

for i=1:8
    arx_model = (B_(1, 1:i+1) * sys_id{i}.B') / (A_(1, 1:i+1) * sys_id{i}.A');
    Y = lsim(arx_model, u);
    for k=1:n
        MSE_id(i) = MSE_id(i) + (y_id(k) - Y(k)) ^ 2;
    end
end
MSE_id = MSE_id ./ n;

% validation data set
MSE_val=zeros(8,1);

for i=1:8
    arx_model = (B_(1, 1:i+1) * sys_val{i}.B') / (A_(1, 1:i+1) * sys_val{i}.A');
    Y_val = lsim(arx_model, u_val);
    for k=1:n
        MSE_val(i) = MSE_val(i) + (y_val(k) - Y_val(k)) ^ 2;
    end
end
MSE_val = MSE_val ./ n;

plot(1:8, MSE_id, 'o-')
hold on;
plot(1:8, MSE_val, 'o-')
title('Mean square error')
xlabel('ARX order')
ylabel('Mean square error')
hold off;
legend('MSE on identification data set', 'MSE on validation data set')

%% compute and plot residuals for different ARX orders
% n = 2,4,6
r_val = zeros(n,8);

for i=1:8
    arx_model = (B_(1, 1:i+1) * sys_val{i}.B') / (A_(1, 1:i+1) * sys_val{i}.A');
    Y_val = lsim(arx_model, u_val);
    for j=1:n
        r_val(j,i) = r_val(j,i) + (y_val(j) - Y_val(j));
    end
end

figure(4)
plot(1:n, r_val(:,2), 'o')
title('Residual n=2')
xlabel('Timestep')
ylabel('Residual')

figure(5)
plot(1:n, r_val(:,4), 'o')
title('Residual n=4')
xlabel('Timestep')
ylabel('Residual')

figure(6)
plot(1:n, r_val(:,6), 'o')
title('Residual n=6')
xlabel('Timestep')
ylabel('Residual')

%% compare frequencies responses
figure(7)
bode(B/A)
hold on;

% delay vectors
A_ = [1 z^-1 z^-2 z^-3 z^-4 z^-5 z^-6 z^-7 z^-8];
B_ = [0 z^-1 z^-2 z^-3 z^-4 z^-5 z^-6 z^-7 z^-8];
for i=[2,4,6]
    bode((B_(1, 1:i+1) * sys_id{i}.B') / (A_(1, 1:i+1) * sys_id{i}.A'))
end

hold off;
legend('True Plant', 'ARX 2', 'ARX 4', 'ARX 6')