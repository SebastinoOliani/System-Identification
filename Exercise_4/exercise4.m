clear all
clc
%%
T=1; %sample time
dim=200;
omega=(2*pi/dim)*[0:dim/2]';
t=T*[0:dim-1]';
% random equally distribuited 1,-1
u1=idinput(dim,'RBS');
% uniform distribuited
u2=rand(dim,1);
U1=fft(u1);
U2=fft(u2);
% noise
v=randn(dim,1)*sqrt(0.01);
%%
z=tf('z',T);
G=(0.0565*z^-1+0.1013*z^-2)/(1-1.4183*z^-1+1.5894*z^-2-1.3161*z^-3+0.8864*z^-4);

y1=lsim(G,u1,t)+v;
y2=lsim(G,u2,t)+v;
Y1=fft(y1);
Y2=fft(y2);

%%
figure(1)
plot(t,y1); hold on; 
plot(t,u1);
title('Equally distribuited +1,-1')
xlabel('time [s]')
ylabel('Magnitude')
legend('y1','u1')
hold off;

figure(2)
plot(t,y2);
hold on;
plot(t,u2);
title('Uniform distribuited inputs')
xlabel('time [s]')
ylabel('Magnitude')
legend('y2','u2')
hold off;

%%
G1=Y1(1:dim/2+1)./U1(1:dim/2+1);
G2=Y2(1:dim/2+1)./U2(1:dim/2+1);
G_freqresponse=squeeze(freqresp(G,omega));

figure(3)
loglog(omega,abs(G_freqresponse))
hold on
loglog(omega,abs(G1))
loglog(omega,abs(G2))
title('ETFE')
xlabel('Omega [rad/s]')
ylabel('Magnitude [dB]')
legend('G(z)','G1(z), equally distr','G2(z), uni distr')
xlim([.2,max(omega)])
ylim([1e-2,1e2])
hold off

figure(4)
loglog(omega,abs(abs(G_freqresponse)-abs(G1)))
hold on
loglog(omega,abs(abs(G_freqresponse)-abs(G2)))
title('Magnitude of the errors')
xlabel('Omega [rad/s]')
ylabel('Magnitude [dB]')
legend('G(z)-G1(z)','G(z)-G2(z)')
xlim([.2,max(omega)])
ylim([1e-4,1e2])
hold off
% random binary signal give smaller error
%%
% longer input
N=1000;
omega=(2*pi/N)*[0:N/2]';
t=T*[0:N-1]';
% random equally distribuited 1,-1
u3=idinput(N,'RBS');
% uniform distribuited
u4=idinput(N,'RBS');
U3=fft(u3);
U4=fft(u4);

y3=lsim(G,u3,t)+randn(N,1)*sqrt(0.01);
y4=lsim(G,u4,t)+mean(randn(N,1)*sqrt(0.01),2);
Y3=fft(y3);
Y4=fft(y4);

G3=Y3(1:N/2+1)./U3(1:N/2+1);
G4=Y4(1:N/2+1)./U4(1:N/2+1);
G_freqresponse=squeeze(freqresp(G,omega));

%repeat esperiment 5 times
u5=idinput(N,'RBS');
y5=lsim(G,repmat(u5,5,1))+randn(N*5,1)*sqrt(0.01);
% discard the first period and average the last four period
y5=mean(reshape(y5(N+1:end),N,[]),2);

U5=fft(u5);
Y5=fft(y5);
G5=Y5(1:N/2+1)./U5(1:N/2+1);

figure(5)
loglog(omega,abs(G_freqresponse))
hold on
loglog(omega,abs(G3))
loglog(omega,abs(G4))
loglog(omega,abs(G5))
title('ETFE')
xlabel('Omega [rad/s]')
ylabel('Magnitude [dB]')
legend('True plant','ETFE 1','ETFE 2','ETFE 3')
xlim([.2,max(omega)])
ylim([1e-2,1e2])
hold off

figure(6)
loglog(omega,abs(abs(G_freqresponse)-abs(G3)))
hold on
loglog(omega,abs(abs(G_freqresponse)-abs(G4)))
loglog(omega,abs(abs(G_freqresponse)-abs(G5)))
title('Magnitude of the errors')
xlabel('Omega [rad/s]')
ylabel('Magnitude [dB]')
legend('G(z)-G3(z)','G(z)-G4(z)','G(z)-G5(z)')
xlim([.2,max(omega)])
ylim([1e-4,1e2])
hold off