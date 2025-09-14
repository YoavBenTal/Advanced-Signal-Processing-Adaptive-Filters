fs = 48000; %Hz
sample_time = 10; %sec

num_samples = sample_time * fs;

%% Question 4
% Generate signal:
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
r = [(1/(1-alpha^2))+sigma_N^2;
    (alpha/(1-alpha^2))]; % this is the first column of the R matrix.

R = toeplitz(r); 

P = [(alpha/(1-alpha^2)); 
    (alpha^2/(1-alpha^2))];

w_opt = R^-1 * P;
norm_W_opt_square = norm(w_opt)^2;

L = 2; % filter order
lambda = 0.99; % forgetting factor, choosing lambda=1 we return to the LMS solution.

%% %% Question 4.a
% RLS algorithm:


Delta = [0.000001, 0.0001 0.1, 10, 1000, 100000];
NR = [];
figure
hold on


for delta = Delta
    Wn = zeros(L,1); % initial guess
    Pn = delta^-1 * eye(L); % initialize Pn matrix.
    X = zeros(num_samples - L,1); % for the plot
    Y = zeros(num_samples - L,1);
    En_vec = zeros(num_samples - L,1);
    for n = L+1:1:num_samples
        % ----- RLS Algo ----------
        Yn = flip(Zn(n-L:n-1));
        Xn_est = transpose(Wn) * Yn; % compute the estimate.

        Xn = Zn(n);
        En = Xn - Xn_est; % compute the error.
        En_vec(n-L,1) = En;

        Kn = (lambda^-1 * Pn * Yn) / (1 + lambda^-1 * transpose(Yn) * Pn * Yn); % compute K vector.

        Wn = Wn + Kn * En; % update Wn.
        Pn = lambda^-1 * Pn - lambda^-1 * Kn * transpose(Yn) * Pn; % update Pn.
        % ----- end ----------


        CnSquare = norm(Wn - w_opt)^2;

        % Plot:

        X(n-L,1) = n;
        Y(n-L,1) = 10*log10(CnSquare / norm_W_opt_square);
    end
    
    NR = [NR 10*log10(sum(Zn.^2)/sum(En_vec.^2))];
    plot(X,Y,"LineWidth",1)
end

title("RLS Algorithm - Normalized error Vs iterations")
xlabel("iterations")
ylabel("Normalized Error dB")
legend("\delta = " + transpose(Delta) + " NR_d_B = " + transpose(NR))
grid on
hold off

