clear all
clc
%%
dim=1024;
e=randn(dim,1);
E_N=fft(e);
per_e=1/dim*abs(E_N).^2;
omega=(2*pi/dim)*[0:dim-1]';
idx = find(omega > 0 & omega < pi); % positive frequencies
figure (1)
loglog(omega(idx),per_e(idx))
title('Periodogram of e')
xlabel('Omega [rad/s]')
ylabel('Periodogram')
%%
%sample time
T=1;
z=tf('z',T);
P=(z+0.5)/(0.5+(z+0.5)*(z-0.5)^2);
t=T*[0:dim-1]';
y=lsim(P,e,t);
figure(2)
plot(t,y)
title('Responde of P(z) to input e(k)')
xlabel('k')
ylabel('y(k)')
%%
Y_N=fft(y);
per_y=1/dim*abs(Y_N).^2;
figure (3)
loglog(omega(idx),per_y(idx))
title('Periodogram of y')
xlabel('Omega [rad/s]')
ylabel('Periodogram')
%%
[mag,phase]=bode(P,omega);
mag=squeeze(mag);
figure(4)
loglog(omega(idx),per_y(idx))
hold on
loglog(omega(idx),abs(mag(idx)).^2)
loglog(omega(idx),abs(per_y(idx)-abs(mag(idx)).^2))
xlabel('Omega [rad/s]')
ylabel('Magnitude [dB]')
legend('1/N*|Y_N(exp(j\omega)|^2','|P(z)|^2','1/N*|Y_N(exp(j\omega)|^2-|P(z)|^2')
hold off
