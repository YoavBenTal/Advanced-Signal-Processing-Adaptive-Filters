fs = 48000; %Hz
sample_time = 10; %sec

num_samples = sample_time * fs;

% Question 1.4

%signal parameters:
sigma_N = 1;
alpha = 0.5;

% Want to make wss signal Zn = Z^n + Nn
Nn = randn(num_samples,1)*sigma_N;

Gn = randn(num_samples,1);
a = [1 -alpha];
b = [1];
Z_tilda_n = filter(b,a,Gn);

% Desired Signal
Zn = Z_tilda_n + Nn;

%% Question 1.4.a
disp("Empirical mean = " + mean(Zn))
disp("Empirical Second Moment = " + (1/num_samples)*sum(Zn.^2))

%% Question 1.4.c
beta = sqrt(3/14);
scaled_Zn = beta*Zn;
% sound(scaled_Zn, fs);
% The amplitude of the signal corresponds to the volume with which to play
% the signal. Hence the need to choose a normalizing factor for volume
% output. In "sound" function it is by default between 1 and -1. We need to
% scale our variance to this in order to avoid clipping if the variance is
% too large, or, alternatively, to avoid unnecessary low volumes if the variance is too low.
% If our signal had non-zero mean value then we would have to correct that
% aswell.

%% Question 1.5.a
% Generate signal with new parameters:
alpha = 0.9;
sigma_N = sqrt(0.5);
Nn = randn(num_samples,1)*sigma_N;
Gn = randn(num_samples,1);
a = [1 -alpha];
b = [1];
Z_tilda_n = filter(b,a,Gn);

% Desired Signal
Zn = Z_tilda_n + Nn;

% calculating R matrix:
r = [(1/(1-alpha^2))+sigma_N^2 (alpha/(1-alpha^2)) (alpha^2/(1-alpha^2)) (alpha^3/(1-alpha^2)) (alpha^4/(1-alpha^2))]; % this is the first row of the R matrix.
p = [(alpha/(1-alpha^2)) (alpha^2/(1-alpha^2)) (alpha^3/(1-alpha^2)) (alpha^4/(1-alpha^2)) (alpha^5/(1-alpha^2))];
W = [];
for L = 1:1:5
    R = toeplitz(r(1:L));
    P = p(1:L)';
    w_n = transpose(R^-1 * P); % these are the optimal linear filters.
    disp("Coefficients for optimal filter of order " + L); disp(""); disp(w_n);
    
    w_n(L+1:5) = 0; % make the filters the same length to organize in matrix.
    W = [W; w_n];
end

%% Question 1.5.b
% Estimate Zn:
Z_estimate_mat = [];
e_n_mat = [];
e_n_scaled_mat = [];
for L = 1:1:5
    b = [0 W(L,1:L)];
    a = [1];
    
    % Z estimator:
    Z_estimator = filter(b,a,Zn);
    Z_estimate_mat = [Z_estimate_mat, Z_estimator];
    
    % Error = Zn - Z estimator
    e_n = Zn - Z_estimator;
    e_n_mat = [e_n_mat, e_n];
    
    % Scale the cancelled signal: en = Zn - Z estimator
    beta = sqrt(1/(2*var(e_n)));
    e_n_scaled = beta.*e_n;
    e_n_scaled_mat = [e_n_scaled_mat, e_n_scaled];
end

%% Question 1.5.c
% Estimation error:
for L = 1:1:5
    mse = (1/num_samples) * sum(e_n_mat(:,L).^2);
    disp("MSE for Zn - Z estimator for L = " + L); disp(""); disp(mse);
end

%% Question 1.5.d

for L = 1:1:5
    
    NR_db = 10*log10(sum(Zn.^2) / sum(e_n_mat(:,L).^2));
    disp("Noise Reduction in dB for L = " + L); disp(""); disp(NR_db);
end



