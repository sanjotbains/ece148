%% HW 8: Lowpass Analog Filter

%% INTRODUCTION
%{
    The objective of this homework assgnment is to contstruct a simple program
    for the design of analog lowpass Butterworth filters.

    The input to my program will be the standard design specifications, including
        (a) passband frequency,
        (b) maximum passband attentuation,
        (c) stopband frequncy, and
        (d) minimum stopband attenuation.

    To test the program, we apply specifications:
        1. Passband frequency:          \omega_p = 5k rad/s
        2. Max. passband attenuation:   a_{max} = 0.5 dB
        3. Stopband frequency:          \omega_s = 10k rad/s
        4. Min. stopband attenuation:   a_{min} = 20 dB

    For the design, my results should include:
        a. the order of the filter,
        b. cutoff frequency omega_c,
        c. locations of the poles,
        d. the transfer function,
        e. frequency-response plot (amplitude only) of the lowpass filter, and
        f. verification of my design by evaluating the frequency response at
            the passband and stopband frequencies
%}      


%% DESIGN SPECIFICATIONS

% Passband Frequency:
omega_p = 5000;     % rad/s
% Maximum Passband Attenuation:
a_max = 0.5;        % dB
% Stopband Frequency:
omega_s = 10_000;   % rad/s
% Minimum Stopband Attenuation:
a_min = 20;         % dB


%% "Step 1: Determine the order \textit{n} of the filter"
%{
    n \geq \frac{ \text{log} \frac{10^{a_{min}/10} - 1}{10^{a_{max}/10} - 1} }
                { 2 \text{log} \frac{\omega_s}{\omega_p} }
%}

numerator = log10(
        (10^(a_min / 10) - 1) ...
      / (10^(a_max / 10) - 1)
    );
denominator = 2 * log10(omega_s / omega_p);

n = ceil(numerator / denominator); % Filter Order
    %   = 5


%% "Step 2: Select the cutoff frequency \omega_0 of the filter"
%{
    \frac { \omega_p } { (10^(a_max / 10) - 1)^{ \frac{1}{2n} } }
    \leq \omega_0 \leq
    \frac { \omega_s } { (10^(a_min / 10) - 1)^{ \frac{1}{2n} } }
%}

min_omega_c =   omega_p / (
                (10^(a_max / 10) - 1)^(1 / (2*n))
            ); %    = 6170.6 rad/s
max_omega_c =   omega_s / (
                (10^(a_min / 10) - 1)^(1 / (2*n))
            ); %    = 6315.9 rad/s

omega_c = 6200;    % Chosen Cutoff Frequency


%% "Step 3: Identify the poles of the transfer function"
%{
    If n is odd:    s = \omega_0 exp( j 2k \pi / 2n )
    If n is even:   s = \omega_0 exp( j (2k+1) \pi / 2n )

    for k = 0, 1, 2, … 2n-1. 
    Select the n poles on the left-half plane from the set of 2n.
%}

% From Step 1, we know n = 5 is odd.
k = 0:(2*n - 1);
s = omega_c * exp( j * k * pi / n );

% The left-half plane is distinguished by the real part of s being negative
s_left = s(real(s) < 0);
%      = -1915.9+5896.6i  -5015.9+3644.3i  -6200.0+0.0i  -5015.9-3644.3i  -1915.9-5896.6i

figure (1);
plot(real(s_left), imag(s_left), 'x', 'LineWidth', 2);
hold on
plot(omega_c*cos(0:0.1:(2*pi)), omega_c*sin(0:0.1:(2*pi)));
axis equal
title('Butterworth Filter Poles in s-plane');
xlabel('Real(s)');
ylabel('Imag(s)');
print('filter_poles.png', '-dpng', '-tight', '-r300');


%% "Step 4: Formulate the transfer function"
%{
    H_n(s) = \frac { \omega_0^n }
                   { (s - s_1) (s - s_2) ... (s - s_n) }
%}

% s = j*omega
H_n = @(omega) omega_c^n ./ prod(bsxfun(@minus, j*omega(:), s_left), 2);


%% Frequency Response Plot (Amplitude)

freq_range = linspace(1, 25000); %  rad/s
freq_response = H_n(freq_range); 

freq_response_dB = 20*log10(abs(freq_response));

% Octave does not have x/yline
xline = @(xval, varargin) line([xval xval], ylim, varargin{:});
yline = @(yval, varargin) line(xlim, [yval yval], varargin{:});
% Usage: xline( 10, 'linestyle', '--', 'color', r' )

figure (2);
plot(freq_range, abs(freq_response), 'color', '#0055ff', 'LineWidth', 2);
hold on
grid on
xline(omega_p, ...
        'LineStyle', '--', ...
        'color', '#747474', ...
        'LineWidth', 1);
text(omega_p, max(ylim)*0.05, '\omega_p ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', '#747474');
xline(omega_c, ...
        'LineStyle', '--', ...
        'color', 'r', ...
        'LineWidth', 1);
text(omega_c, max(ylim)*0.05, '\omega_c ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', 'r');
xline(omega_s, ...
        'LineStyle', '--', ...
        'color', '#747474', ...
        'LineWidth', 1);
text(omega_s, max(ylim)*0.05, '\omega_s ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', '#747474');
title('Butterworth Filter Frequency Response');
ylabel('Magnitude Gain (scalar)');
xlabel('Frequency (rad/s)');
print('freq_response.png', '-dpng', '-tight', '-r300');

figure (3);
semilogx(freq_range, freq_response_dB, 'color', '#0055ff', 'LineWidth', 1.5);
ylim([-25 2.5]);
hold on
xline(omega_p, 'LineStyle', '--','color', '#747474', 'LineWidth', 1);
text(omega_p, min(ylim)*0.95, '\omega_p ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', '#747474');
yline(-0.5, 'LineStyle', '--','color', '#747474', 'LineWidth', 1)
xline(omega_s, 'LineStyle', '--', 'color', '#747474', 'LineWidth', 1);
text(omega_s, min(ylim)*0.95, '\omega_s ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', '#747474');
yline(-20, 'LineStyle', '--','color', '#747474', 'LineWidth', 1)
grid on
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Butterworth Filter Frequency Response (dB)');
print('filter_bode_plot.png', '-dpng', '-tight', '-r300');


% The attenuation at the passband frequency is:
omega_p_dB = 20*log10(abs(H_n(omega_p))); %     = -0.4780 dB
% The attenuation at the stopband frequency is:
omega_s_dB = 20*log10(abs(H_n(omega_s))); %     = -20.797 dB
