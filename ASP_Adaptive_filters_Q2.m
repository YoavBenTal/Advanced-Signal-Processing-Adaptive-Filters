fs = 48000; %Hz
sample_time = 10; %sec

num_samples = sample_time * fs;

%% Question 2.a
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

L = 4; % filter length

% calculating R matrix:
r = [(1/(1-alpha^2))+sigma_N^2;
    (alpha/(1-alpha^2));
    (alpha^2/(1-alpha^2));
    (alpha^3/(1-alpha^2))]; % this is the first column of the R matrix.

R = toeplitz(r); 

P = [(alpha/(1-alpha^2)); 
    (alpha^2/(1-alpha^2));
    (alpha^3/(1-alpha^2));
    (alpha^4/(1-alpha^2))];

w_opt = R^-1 * P;
eigen_values = eigs(R)
% The eigen values are important: we need -1 < 1-mu*lambda < 1 where lambda
% is the largest eigenvalue for this algorithm to converge.

%% Question 2.b
% Steepest Descent Algorithm and error calculations: 
figure
hold on

Mu = [0.001 0.01 0.1 0.2];

for mu = Mu
    Wn = [0;0;0;0]; % initial guess
    
    X = []; % for the plot
    Y = [];
    for n = 0:1:99  % 100 iterations of steepest descent
        Wn = Wn + mu * (P - R*Wn);
        CnSquare = norm(Wn - w_opt)^2;
        
        % Plot:
        
        X = [X; n];
        Y = [Y; 10*log10(CnSquare / norm(w_opt)^2)];
    end
    
    plot(X,Y,"LineWidth",2)
    
end

title("Steepest Descent Algo - Normalized error Vs iterations")
xlabel("iterations")
ylabel("Normalized Error")
legend("\mu = " + Mu)
ylim([-100 30])
grid on
hold off
