set(0,'DefaultFigureVisible','off'); % Set default visibility to off
% Set default font to match LaTeX serif font
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

%% Background

% We have a periodic signal f(t) with period T, w/ 7 harmonics of the form:
%   f(t) = \sum_{m=1}^{7} a_m * sin(m * w_0 * t)
% where w_0 = 2 * pi / T is the fundamental frequency.

T = 1; % Period of the signal
w_0 = 2 * pi / T; % Fundamental frequency
a_m = [9 5 8 7 1 8 9]; % Amplitude of harmonics

f = @(t) a_m(1) * sin(1 * w_0 * t) + a_m(2) * sin(2 * w_0 * t) + ...
         a_m(3) * sin(3 * w_0 * t) + a_m(4) * sin(4 * w_0 * t) + ...
         a_m(5) * sin(5 * w_0 * t) + a_m(6) * sin(6 * w_0 * t) + ...
         a_m(7) * sin(7 * w_0 * t);

% We will plot the function f(t) over one period.
t = linspace(0, T, 1000); % Time vector for one period
f_t = f(t); % Function values over one period

% We then take 32 uniform samples of the function f(t) over one period:
delta_T = T / 32; % Sampling interval
t_samples = 0:delta_T:(T-delta_T); % Sampled time vector
f_n = f(t_samples); % Sampled function values {f(n = 0, 1, ..., 31)}

% To observe the spectrum of the function f(t), we will take a 32-point DFT 
% of the sampled function values. 
F_k = fftshift(fft(f_n)); % F(k) = DFT of f(n = 0, 1, ..., 31)

%% Problem 1: Time Function and Sampling
%{
    Plot one full period of the function f(t), over 0 <= t < T
     and the sampled values f(n = 0, 1, ..., 31).
%}

% Plot the periodic signal f(t) over one period
f_t_plot = figure (1); % Create a new figure
    plot(t, f_t, 'r', 'LineWidth', 2);
    hold on;
    stem(t_samples, f_n, 'b', 'LineWidth', 2); 
    xlabel('Time (t)');
    ylabel('f(t)');
    grid on;
    axis([0 T -36 36]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(f_t_plot, 'f_t.png', 'Resolution', 300);

%% Problem 2: Fourier Series Expansion
%{
    Formulate the periodic signal f(t) in the form of complex Fourier series expansion. 
    List and sketch the Fourier coefficients {F_m} as a sequence of m. 
    And then identify the corresponding physical frequencies.
%}

N = 7;  % Number of harmonics to compute (from -N to N)
m = -N:N;  % Harmonic indices corresponding to F_m 
F_m = zeros(size(m));  % Complex Fourier coefficients (vector from -N to N) 
dt = t(2) - t(1);   % Calculate dt for numerical integration
    
% Calculate each coefficient using rectangular integration
for i = 1:length(m)
    integrand = f_t .* exp(-1j * m(i) * w_0 * t);
    F_m(i) = (1/T) * sum(integrand) * dt; 
end
clear i integrand;

fprintf('Fourier Coefficients:\n');
for i = 1:length(m)
    fprintf('F_%d = %.1f + %.1fj\n', m(i), real(F_m(i)), imag(F_m(i)));
end
clear i;

% Plot the Fourier coefficients
F_m_plot = figure (2); % Create a new figure
    stem(m, abs(F_m), 'LineWidth', 2);
    xlabel('Harmonic Index (m)');
    ylabel('Magnitude of Fourier Coefficients');
    grid on;
    axis([-N-1 N+1 -0.5 5]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(F_m_plot, 'F_m.png', 'Resolution', 300);

% Calculate the physical frequencies corresponding to the Fourier coefficients
f_m = m * w_0; % Physical frequencies
fprintf('Physical Frequencies:\n');
for i = 1:length(m)
    fprintf('f_%d = %d Hz\n', m(i), f_m(i)/(2*pi));
end

%% Problem 3: Fourier Transform of Periodic Functions
%{
    Determine and sketch the Fourier transform F(jω) of the periodic function f(t). 
    Identify the amplitudes and physical frequencies of the peaks.
%}

% sin(ωt) = (e^(jωt) - e^(-jωt)) / (2j)
% ω = [1:7] * w_0
% F(jω) = 1/2j * (F_1 * δ(ω - ω_1) - F_-1 * δ(ω + ω_1) + 
%         ... + F_7 * δ(ω - ω_7) - F_-7 * δ(ω + ω_7))

% Plot the Fourier transform F(jω)
F_j_omega_plot = figure (3); % Create a new figure
    stem(f_m, abs(F_m), 'LineWidth', 2);
    xlabel('Frequency (ω)');
    ylabel('Magnitude of Fourier Transform');
    xticks([-N*2*pi:2*pi:N*2*pi]); % Set x-axis ticks in units of pi
    xticklabels({'-14π', '-12π', '-10π', '-8π', '-6π', '-4π', '-2π', '0', '2π', '4π', '6π', '8π', '10π', '12π', '14π'});
    grid on;
    axis([-N*2*pi N*2*pi -0.5 5]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(F_j_omega_plot, 'F_j_omega.png', 'Resolution', 300);

%% Problem 4: DTFT
%{
    Plot the DTFT F(e^jΩ) of the 32-point time-domain sequence f(n) over the interval (-π, +π). 
    Identify the amplitudes and locations of the peaks
%}

% Define the frequency range (-π to π)
Omega = linspace(-pi, pi, 1000);

% Initialize the DTFT result
F_Omega = zeros(size(Omega));
% Compute the DTFT using the Fourier coefficients
for n = 1:length(Omega)
    F_Omega(n) = sum(f_n .* exp(-1j * Omega(n) * (0:31)));
end

% Plot the DTFT result
F_Omega_plot = figure (4); % Create a new figure
    plot(Omega, real(F_Omega), 'LineWidth', 2);
    xlabel('Frequency (Ω)');
    ylabel('Magnitude of DTFT');
    grid on;
    axis([-pi pi -150 230]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(F_Omega_plot, 'F_Omega.png', 'Resolution', 300);

%% Problem 5: DFT
%{ 
    Compute and plot the 32-point DFT of this 15-point spectral sequence F(k).
    Identify the amplitudes and locations of the peaks.
%}

% We computed the DFT in the background section
% F(k) = DFT of f(n = 0, 1, ..., 31)

% Plot the DFT result
F_k_plot = figure (5); % Create a new figure
    plot(-16:15, real(F_k), 'LineWidth', 2);
    xlabel('Frequency Bin (k)');
    ylabel('Magnitude of DFT');
    grid on;
    axis([-16 15 -inf inf]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(F_k_plot, 'F_k.png', 'Resolution', 300);
