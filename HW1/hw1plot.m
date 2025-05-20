% Parameters
N = 5;                        % Number of positive/negative harmonics
omega_x = 1;                 % Unit frequency (can be symbolic or actual)
n = -N:N;                    % Indices
frequencies = n * omega_x;   % Impulse locations
amplitudes = rand(1, 2*N+1); % Random amplitudes

%%
% Plot
figure (1);
stem(frequencies, amplitudes, 'filled', 'LineWidth', 1.2, 'Color', 'r');
xlabel('\omega (in units of \omega_x)');
ylabel('Amplitude');
title(['Fourier Spectrum of g(t) for N = ', num2str(N)]);
xticks(frequencies);
xticklabels(arrayfun(@(x) sprintf('%d\\omega_x', x), n, 'UniformOutput', false));
xlim([-7 7]);
set(gca,'ytick',[])
