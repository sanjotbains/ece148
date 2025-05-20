% We have a periodic signal f(t) with period T = 16:
f = @(t) - 5 * (mod(t+8, 16)-8 >= -8  &  mod(t+8, 16)-8 < -5) ...
         + 3 * (mod(t+8, 16)-8 >= -5  &  mod(t+8, 16)-8 <  5) ...
         - 5 * (mod(t+8, 16)-8 >=  5  &  mod(t+8, 16)-8 <  8);

T = 16; % Period of the signal
t = linspace(-T/2, T/2, 1000); % Time vector over one period

% Set default font to match LaTeX serif font
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

%% Problem 1
%{
    Plot the function f(t) and formulate the Fourier series expansion of the function, in the complex form.
%}

% Plot the periodic signal f(t) over one period
f_t_plot = figure (1); % Create a new figure
    plot(t, f(t), 'r', 'LineWidth', 2);
    xlabel('Time (t)');
    ylabel('f(t)');
    grid on;
    axis([-T/2 T/2 -7 5]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(f_t_plot, 'f_t.png', 'Resolution', 300);

% Fundamental frequency
w_0 = 2*pi/T;

% Number of harmonics to compute (from -N to N)
N = 7; 

% Harmonic indices corresponding to c_k
k = -N:N;   

% Complex Fourier coefficients (vector from -N to N)
F_n = zeros(size(k));   
    
% Calculate dt for numerical integration
dt = t(2) - t(1);
    
% Compute the function values
f_t = f(t);
    
% Calculate each coefficient using rectangular integration
for n = 1:length(k)
    integrand = f_t .* exp(-1j * k(n) * w_0 * t);
    F_n(n) = (1/T) * sum(integrand) * dt; 
end

%% Problem 2
%{
    For this periodic signal, utilize 7 harmonics, for n = ± 1, ± 2, ….± 7, to reconstruct the signal. Plot the result for one period, within the interval (-8, +8).
%}

% Reconstruct the function from Fourier series
f_reconstructed = zeros(size(t));
for n = 1:length(k)
    f_reconstructed = f_reconstructed + F_n(n) * exp(1j * k(n) * w_0 * t);
end

% Plot the reconstructed function
f_rec_plot = figure (2); % Create a new figure
    plot(t, real(f_reconstructed), 'LineWidth', 2);
    % title('Reconstructed Function from Fourier Series');
    hold on; % Hold the current plot
    plot(t, f(t), 'r--', 'LineWidth', 2); % Plot original function for comparison
    legend('Reconstructed', 'Original', 'Location', 'south'); % Add legend
    hold off; % Release the hold on the current plot
    xlabel('Time (t)');
    ylabel('f(t)');
    grid on;
    axis([-T/2 T/2 -7 5]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(f_rec_plot, 'f_reconstructed.png', 'Resolution', 300);

%% Problem 3
%{
    The 15 Fourier coefficients corresponding to the 7 harmonics are grouped into a sequence {Fn}, where n = -7, -6, … +5, +6, +7. 
    Take the DTFT of this 15-point sequence and plot the result within the interval (-π, +π).
%}

% Define the frequency range (-π to π)
omega = linspace(-pi, pi, 1000);

% Initialize the DTFT result
F_omega = zeros(size(omega));
% Compute the DTFT using the Fourier coefficients
for k = 1:length(omega)
    F_omega(k) = sum(F_n .* exp(-1j * omega(k) * (-N:N)));
end

% Plot the DTFT result
DTFT_plot = figure (3); % Create a new figure
    plot(omega, real(F_omega), 'LineWidth', 2); % Plot the magnitude of the DTFT
    % title('DTFT of Fourier Coefficients');
    xlabel('Frequency (\omega)');
    ylabel('Re\{F(\omega)\}');
    xticks([-pi -pi/2 0 pi/2 pi]); % Set x-axis ticks in units of pi
    xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'}); % Set corresponding labels
    grid on;
    axis([-pi pi -7 5]);
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(DTFT_plot, 'DTFT.png', 'Resolution', 300); 

%% Hints:
%{
    For Problems 4 and 5 of Assignment 4, the sequence you use as the input to the FFT should be in the form of
    {0, F_1, F_2, F_3, F_4, F_5, F_6, F_7, 0, ..., 0, F_-7, F_-6, F_-5, F_-4, F_-3, F_-2, F_-1}
    You need to insert 17 zeros to the middle portion to make it 32 points for Problem 4 and 49 zeros for Problem 5 to make it 64 points.
    The first term, which is zero, is the DC term, since the DC level of the periodic function is zero. 
    The rearrangement of the Fourier coefficients into this form is because of the index system of FFT. 
    When you complete the FFT, you will need to reverse the arrangement to visualize the time function.
%}

%% Problem 4
%{
    Take the 32-point DFT of this 15-point sequence {F_n} and plot the result.
%}

% Create a 32-point sequence with zeros in the middle
F_32 = zeros(1, 32); % Initialize a 32-point sequence with zeros
F_32(2:8) = F_n(9:15); % Fill in the positive Fourier coefficients (F_1 to F_7)
F_32(26:32) = F_n(1:7); % Fill in the negative coefficients (F_-7 to F_-1)

F_fft_32 = fft(F_32); % Compute the 32-point DFT
F_fft_32_shifted = fftshift(F_fft_32); % Shift the zero frequency component to the center

% Plot the real part of the DFT result
DFT_32_plot = figure (4); % Create a new figure
    plot(-16:15, real(F_fft_32_shifted), 'LineWidth', 2); % Plot the magnitude of the DFT with binning from -N/2 to N/2
    % title('32-point DFT of Fourier Coefficients');
    xlabel('Frequency Bins (k)');
    ylabel('Real\{F(k)\}');
    grid on;
    axis([-16 15 -7 5]); % Adjust axis limits to match the binning
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(DFT_32_plot, 'DFT_32.png', 'Resolution', 300);

%% Problem 5
%{
    Take the 64-point DFT of this 15-point sequence {F_n} and plot the result.
%}

F_64 = zeros(1, 64); % Initialize a 64-point sequence with zeros
F_64(2:8) = F_n(9:15); % Fill in the positive Fourier coefficients (F_1 to F_7)
F_64(58:64) = F_n(1:7); % Fill in the negative coefficients (F_-7 to F_-1)

F_fft_64 = fft(F_64); % Compute the 64-point DFT
F_fft_64_shifted = fftshift(F_fft_64); % Shift the zero frequency component to the center

% Plot the real part of the DFT result
DFT_64_plot = figure (5); % Create a new figure
    plot(-32:31, real(F_fft_64_shifted), 'LineWidth', 2); % Plot the magnitude of the DFT with binning from -N/2 to N/2
    % title('64-point DFT of Fourier Coefficients');
    xlabel('Frequency Bins (k)');
    ylabel('Real\{F(k)\}');
    grid on;
    axis([-32 31 -7 5]); % Adjust axis limits to match the binning
    set(gca, 'FontName', 'Times New Roman'); % Apply font to axes
exportgraphics(DFT_64_plot, 'DFT_64.png', 'Resolution', 300);
