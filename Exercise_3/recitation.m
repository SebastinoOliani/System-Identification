clear all
clc

%%
N=400;
t=0:N-1;
u=1.*randn(N,1);
figure (1)
plot(t,u,'.-')
title('Input signal u')
xlabel('k [sample]')
ylabel('u')
legend('u')

% autocorrelation of u
% positive time lags
lags=0:(N-1)/2;
%lags=-N/2+1:N/2;
R=zeros(N,1);
R(1,1)=u(1)*u(1)/N;
for i=1:1:200
    for k=2:1:N
        R(k,1)=u(k)*u(k)/N;
    end
end
%R=autocorr(u(1:(N)/2,1),N/2-1);
omega=-N/2+1:N/2;
figure(2)
plot(omega,R,'.-')
title('Periodic autocorrelation R')
xlabel('lags [Sample]')
ylabel('R')
xlim([-N/2,N/2])
ylim([-0.2,1])
legend('R')

%%
N=400;
t=0:N-1;
u=sawtooth(t.*2*pi/50,1/2)';
figure (1)
plot(t,u,'.-')
title('Input signal u')
xlabel('k [sample]')
ylabel('u')
legend('u')

% autocorrelation of u
lags=-N/2+1:N/2;
R=zeros(N,1);
for k=1:1:N
    R(k,1)=u(k)*u(k-lags(1,k))/N;
end
figure(2)
plot(lags,R,'.-')
title('Periodic autocorrelation R')
xlabel('lags [Sample]')
ylabel('R')
legend('R')

%%
N=400;
t=0:N-1;
u=sawtooth(t.*2*pi/50,1/2)'+0.2.*randn(N,1);
figure (1)
plot(t,u,'.-')
title('Input signal u')
xlabel('k [sample]')
ylabel('u')
legend('u')

%%
N=430;
t=0:N-1;
u=sawtooth(t.*2*pi/50,1/2)';
figure (1)
plot(t,u,'.-')
title('Input signal u')
xlabel('k [sample]')
ylabel('u')
legend('u')