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

F_k = fftshift(fft(f_n)); % F(k) = DFT of f(n = 0, 1, ... , 15)

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

F_a = fftshift(fft(f_a));

F_a_plot = figure (1); % Create a new figure
    plot(-32/64:1/64:31/64, abs(F_a), 'LineWidth', 2);
    hold on
    plot(-32/64:4/64:31/64, abs(F_k), 'LineWidth', 2);
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
    48 zeros in the middle. The extended spectral sequence 
%}