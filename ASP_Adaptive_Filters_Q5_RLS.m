%%%% Q5 RLS Algorithm %%%%
[air_plane, Fs] = audioread('airplane.wav');
[city, ] = audioread('city.wav');
[cafe, ] = audioread('cafe.wav');
[vac_clean, ] = audioread('vacuumcleaner.wav');

num_samples = length(city);


L = 4; % filter order
lambda = 0.99; % forgetting factor, choosing lambda=1 we return to the LMS solution.
delta = 0.1;

M = 10000; % num of samples for instantenouse power.

NR = [];

sig = air_plane;
Wn = zeros(L,1); % initial guess
Pn = delta^-1 * eye(L); % initialize Pn matrix.

En_vec = zeros(num_samples,1);
IP_sig_vec = zeros(num_samples - M,1);
IP_filtered_vec = zeros(num_samples - M,1);
    for n = L+1:1:num_samples
        % ----- RLS Algo ----------
        Yn = flip(sig(n-L:n-1));
        Xn_est = transpose(Wn) * Yn; % compute the estimate.

        Xn = sig(n);
        En = Xn - Xn_est; % compute the error.
        En_vec(n-L,1) = En;

        Kn = (lambda^-1 * Pn * Yn) / (1 + lambda^-1 * transpose(Yn) * Pn * Yn); % compute K vector.

        Wn = Wn + Kn * En; % update Wn.
        Pn = lambda^-1 * Pn - lambda^-1 * Kn * transpose(Yn) * Pn; % update Pn.
        % ----- end ----------

        % calc inst power:
        if n > M
            IP_sig = 10*log10((1/M) * sum(sig(n-M+1:n).^2));
            IP_sig_vec(n-M, 1) = IP_sig;
        end
        if n > M+L+1
            IP_filtered = 10*log10((1/M) * sum(En_vec(n-M+1:n).^2));
            IP_filtered_vec(n-M,1) = IP_filtered;
            
        end

    end
    
    NR = [NR 10*log10(sum(sig.^2)/sum(En_vec(30000:end).^2))];
    X = transpose(M:1:num_samples -1);
    
figure
plot(X,IP_sig_vec,"LineWidth",1)
hold on
plot(X,IP_filtered_vec,"LineWidth",1)


title("airplane.wav, RLS, L = " + L + " \lambda = " + lambda + " \delta = " + delta + " NR = " + NR)
xlabel("iterations")
ylabel("Instantenouse Power [dB]")
legend("Original signal", "Filtered Signal")
grid on
hold off

