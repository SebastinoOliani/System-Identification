function [R_u,lags,omegas,U_cos,U_cos_Wf,U_cos_Wt] = HS2022_SysID_Exercise_05_22952683()

%% Output format specification
% R_u must be a 3xN matrix
% lags must be a 1xN vector
% omegas must be a 1xN vector
% U_cos must be a 1xN vector
% U_cos_Wf must be a 1xN vector
% U_cos_Wt must be a 1xN vector
%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

[u_prbs,u_rand,u_cos] = HS2022_SysID_Exercise_05_GenerateData(LegiNumber);

N = length(u_prbs);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the input signals u_prbs, u_randn and u_cos to solve the problem.

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% 1. Calculation of autocorrelation


function [R, lags] = autocorrelation(N, u)

  R = zeros(N,1);
  
  for k = 0:N-1
    u_shifted = [u(k+1:end), u(1:k)];
    R(k+1,1) = 1/N * sum(u .* u_shifted);
  end

  lags = -(N-1)/2:(N-1)/2;

  % Shift values to put zero in the middle
  mid = (N+1)/2;
  R = [R(mid+1:end,1); R(1:mid,1)];

end

R_u = NaN * ones(3,N);
[R1,lags1]=autocorrelation(N,u_prbs);
[R2,lags2]=autocorrelation(N,u_rand);
[R3,lags3]=autocorrelation(N,u_cos);
R_u(1,:)=R1;
R_u(2,:)=R2;
R_u(3,:)=R3;

figure(1)
plot(lags1,R_u(1,:),'.-');
hold on
plot(lags2,R_u(2,:),'.-');
plot(lags3,R_u(3,:),'.-');
title('Periodic autocorrelation of the input signals')
xlabel('lags [Sample]')
ylabel('Autocorrelation')
legend('R_(prbs)','R_(rand)','R_(cos)')
xlim([-N/2+1;N/2])
ylim([-0.5;1])
hold off;


%% 2. Smoothing

omegas = (2*pi/N)*[-(N-1)/2:(N-1)/2]';
%omegas = omegas-pi;
U_cos = fft(u_cos);
U_cos = [U_cos((N+1)/2+1:end),U_cos(1:(N+1)/2)];

figure (2)
plot(omegas,abs(U_cos));
title('Magnitude of U_(cos)')
xlabel('Omega [rad/s]')
ylabel('Magnitude')
legend('U_(cos)')

% Hann windows
gamma=150;
U_cos_Wf = 0 * ones(1,N);

%results are in [-pi;pi]
[omega,Wf_Hann]=WfHann(gamma,N);
z_idx=find(omega==0);

for wn=1:N
    Wnorm=0;
    for xi=1:N
        widx=mod(xi-wn,N)+1;
        U_cos_Wf(wn)=U_cos_Wf(wn)+Wf_Hann(widx)*U_cos(xi);
        Wnorm=Wnorm+Wf_Hann(widx);
    end
    U_cos_Wf(wn)=U_cos_Wf(wn)/Wnorm;
end
U_cos_Wf=[U_cos_Wf(z_idx+1:N),U_cos_Wf(1:z_idx)];

U_cos_Wt = NaN * ones(1,N);
[lags,Wt_Hann]=WtHann(gamma,N);
z_idx=find(lags==0);
u_cos_smoothed=u_cos.*Wt_Hann';
U_cos_Wt=fft(u_cos_smoothed);
U_cos_Wt=[U_cos_Wt(z_idx+1:N),U_cos_Wt(1:z_idx)];
lags=2*pi/N*lags;

figure (3)
plot(lags,abs(U_cos_Wt),'.-','LineWidth',2);
hold on;
plot(omega,abs(U_cos_Wf),'.-');
title('Magnitude of UcosWf and UcosWt')
xlabel('Omega [rad/s]')
ylabel('Magnitude')
legend('Time domain smoothing','Frequency domain smoothing')
hold off;

%% Functions

    function [lags_w,wHann] = WtHann(gamma,N)
        %------------------------------------------------------------------
        %
        %   [lags,wHann] = WtHann(gamma,N)
        %
        %   Create a Hann window with width parameter gamma and data length N.
        %   The Hann window is a raised cosine.
        %
        %   Roy Smith,  18 October, 2017.
        %
        %------------------------------------------------------------------

        if nargin == 0,
            disp('Syntax: [lags,w] = WtHann(gamma,N)')
            return
        elseif nargin ~= 2,
            error('incorrect number of input arguments (2 expected)')
            return
        end

        %   basic parameter checking
        if length(gamma) > 1,
            error('Width parameter, gamma, must be a scalar');
        end
        if round(gamma) ~= gamma,
            error('Width parameter, gamma, must be an integer');
        end
        if gamma < 1,
            error('Width parameter, gamma, must be positive');
        end
        if length(N) > 1,
            error('Calculation length, N, must be a scalar');
        end
        if round(N) ~= N,
            error('Calculation length, N, must be an integer');
        end
        if N < 1,
            error('Calculation length, N, must be positive');
        end

        lags_w = [floor(-N/2+1):floor(N/2)]';
        wHann = 0*lags_w;
        idx = find(abs(lags_w) <= gamma);
        wHann(idx) = 0.5*(1+cos(pi*lags_w(idx)/gamma));

    end

%--------

    function [omega,WHann] = WfHann(gamma,N)
        %------------------------------------------------------------------
        %
        %   [omega,WHann] = WfHann(gamma,N)
        %
        %   Create a frequency domain Hann window with width parameter gamma
        %   and data length N.  The Hann window is a raised cosine.
        %
        %   Roy Smith,  18 October, 2017.
        %
        %                6 November, 2017.  Fixed bug in N even indexing.
        %
        %------------------------------------------------------------------

        if nargin == 0,
            disp('Syntax: [omega,W] = WfHann(gamma,N)')
            return
        elseif nargin ~= 2,
            error('incorrect number of input arguments (2 expected)')
            return
        end

        %   basic parameter checking
        if length(gamma) > 1,
            error('Width parameter, gamma, must be a scalar');
        end
        if round(gamma) ~= gamma,
            error('Width parameter, gamma, must be an integer');
        end
        if gamma < 1,
            error('Width parameter, gamma, must be positive');
        end
        if length(N) > 1,
            error('Calculation length, N, must be a scalar');
        end
        if round(N) ~= N,
            error('Calculation length, N, must be an integer');
        end
        if N < 1,
            error('Calculation length, N, must be positive');
        end

        %   The simplest approach is to define the window in the time domain and
        %   then transform it to the frequency domain.

        lags_w = [floor(-N/2+1):floor(N/2)]';
        wHann = 0*lags_w;
        idx = find(abs(lags_w) <= gamma);
        wHann(idx) = 0.5*(1+cos(pi*lags_w(idx)/gamma));

        %
        zidx = find(lags_w==0);    % index of the zero point.

        wH_raw = fft([wHann(zidx:N);wHann(1:zidx-1)]);
        WHann(zidx:N) = wH_raw(1:N-zidx+1);  % shift +ve freq to end
        WHann(1:zidx-1) = wH_raw(N-zidx+2:N);% shift > pi freq to beginning
        WHann = real(WHann);   % should actually be real
        omega = 2*pi/N*lags_w;

    end

end
