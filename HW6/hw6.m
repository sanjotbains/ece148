set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

%% Interpolation by DFT
%{
    Consider a real periodic signal f(t) of period T. The signal has 7
    harmonics, for m = 1, 2, 3, ... , 7.

        f(t) = \sum_{m=1}^7 a_m \text{sin}(m \omega_0 t)

    where \omega_0 = \frac{2\pi}{T} and the 7 coefficients {a_m} are 
    formed with my 7-digit perm number: 9587189.
%}

T = 1; % Period of the signal
w_0 = 2 * pi / T; % Fundamental frequency
a_m = [9 5 8 7 1 8 9]; % Amplitude of harmonics

f = @(t) a_m(1) * sin(1 * w_0 * t) + a_m(2) * sin(2 * w_0 * t) + ...
         a_m(3) * sin(3 * w_0 * t) + a_m(4) * sin(4 * w_0 * t) + ...
         a_m(5) * sin(5 * w_0 * t) + a_m(6) * sin(6 * w_0 * t) + ...
         a_m(7) * sin(7 * w_0 * t);

%{
    Then we take 16 uniform samples within one period with sample spacing
    \Delta t = \frac{T}{16} to form a short 16-point sequence {f(n)}, 
    where n = 0, 1, 2, ... , 15.
%}

delta_T = T / 16; % Sampling interval
t_samples = 0:delta_T:(T-delta_T); % Sampled time vector
f_n = f(t_samples); % Sampled function values {f(n = 0, 1, ... , 15)}

%{
    Subsequently, we take a 16-point DFT of the sequence to obtain the 
    16-point spectral sequence F(k), where k = 0, 1, 2, ... , 15.

        F(k) = DFT_{N = 16} {f(n)}
%}

F_k = fft(f_n); % F(k) = DFT of f(n = 0, 1, ... , 15)
F_k_shifted = fftshift(F_k);

%% 1. Interpolation of the DFT Spectrum
%{
    Extend the sequence f(n) to 64 points by padding 48 zeros. The 
    extended sequence f_a(n) is in the form

        f_a(n) = f(n)   n = 0, 1, 2, ... , 15
                 0      n = 16, 17, 18, ... , 63
%}

f_a = zeros(1, 64);
f_a(1:16) = f_n;

%{
    Compute and plot the 64-point DFT F_a(k). Compare F_a(k) with F(k) 
    and summarize your observations.

        F_a(k) = DFT_{N = 64} {f_a(n)}
%}

F_a_shifted = fftshift(fft(f_a));

F_a_plot = figure (1); % Create a new figure
    plot(-32/64:1/64:31/64, abs(F_a_shifted), 'b', 'LineWidth', 2);
    hold on
    plot(-32/64:4/64:31/64, abs(F_k_shifted), 'r', 'LineWidth', 2);
    xlabel('Normalized Frequency (\omega/2\pi)');
    ylabel('Magnitude of DFT');
    legend('64-point DFT', '16-point DFT', 'Location', 'best');
    grid on;
    axis([-0.5 0.5 -inf inf]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(F_a_plot, 'F_a.png', 'Resolution', 300);

%% 2. Interpolation in Time Domain
%{
    Extend the 16-point spectral sequence F(k) to 64 points by inserting 
    48 zeros in the middle. The extended spectral sequence is denoted as
    F_b(k).
%}

F_b = zeros(1, 64);
F_b(1:8) = F_k(1:8);
F_b(57:64) = F_k(9:16);
F_b = 4 * F_b;

%{
    Perform a 64-point inverse DFT to bring it back to the time domain

        f_b(n) = IDFT_{N = 64} {F_b(k)}

    Plot the 64-point sequence f_b(n). Compare f_b(n) with f(n) and 
    summarize your observations.
%}

f_b = ifft(F_b);

f_b_plot = figure (2);
    n_64 = 0:0.25:15.75;
    n_16 = 0:15;

    plot(n_64, real(f_b), 'b-', 'LineWidth', 2); hold on;
    stem(n_16, real(f_n), 'r', 'filled', 'LineWidth', 1.5);

    xlabel('n');
    ylabel('Amplitude');
    legend('Interpolated f_b(n)', 'Original f_n', 'Location', 'best');
    grid on;
    set(gca, 'FontName', 'Times New Roman');
exportgraphics(f_b_plot, 'f_b.png', 'Resolution', 300);

%% Signal Scrambling
%{
    The objective of this exercise is to implement a simple digital 
    speech scrambler.

    We use the microphone of of our computer to record a short speech 
    signal g(t), and digitize the speech signal with the A/D tool in
    Audacity into the discrete form g(n).
%}

[g_n, f_s] = audioread('g_n.wav');

%% 1. DFT Spectrum
%{
    Display the DFT spectrum G(k) of the digitized speech signal g(n).
%}

G_k = fft(g_n);
G_k_shifted = fftshift(G_k);

N = length(g_n);
f_axis = linspace(-f_s/2, f_s/2, N);

G_k_plot = figure (3);
    spectrum_db = 20*log10(abs(G_k_shifted) + eps);
    spectrum_smooth = smoothdata(spectrum_db, 'gaussian', 100);
    
    plot(f_axis, spectrum_db, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); 
    hold on;
    plot(f_axis, spectrum_smooth, 'b', 'LineWidth', 1);

    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Magnitude Spectrum of g(n)');
    legend('Raw Spectrum', 'Smoothed Spectrum', 'Location', 'best');
    grid on;

    set(gca, 'XScale', 'log');
    set(gca, 'FontName', 'Times New Roman');
exportgraphics(G_k_plot, 'G_k.png', 'Resolution', 300);

%% 2. Speech Scrambling
%{
    Apply the speech scrambling procedure to the digitized speech signal 
    g(n) and display the DFT spectrum of the scrambled speech signal 
    ĝ(n). Then use the D/A tool to convert it back to an analog signal 
    to check if it is audible.
%}
    
scrambling_seq = ones(size(g_n));
scrambling_seq(2:2:end) = -1;

g_hat_n = g_n .* scrambling_seq;

G_hat_k = fft(g_hat_n);
G_hat_k_shifted = fftshift(G_hat_k);

figure(4);
    spectrum_db = 20*log10(abs(G_hat_k_shifted) + eps);
    spectrum_smooth = smoothdata(spectrum_db, 'gaussian', 100);
    
    plot(f_axis, spectrum_db, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); 
    hold on;
    plot(f_axis, spectrum_smooth, 'b', 'LineWidth', 1);
    
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Spectrum of Scrambled Speech');
    grid on;
    set(gca, 'XScale', 'log');

% Save scrambled audio
audiowrite('g_hat_n.wav', g_hat_n, f_s);

%% 3. Descrambling
%{
    The procedure for scrambling a discrete sequence is simply a 
    multiplication process by the sequence. And we perform the 
    descrambling process with the same sequence. It is common that the 
    scrambling-descrambling process is not exactly synchronized and the 
    offset produces an extra {-1} factor. It results in -g(t), instead 
    of g(t). Check if it is audible when the offset occurs.
%}