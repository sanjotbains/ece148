N = 128;
M = 512;

[m, n] = meshgrid(0:N-1, 0:M-1);  % Grid of indices
f = exp(1j * 2 * pi * m .* n / N);  % Compute f(m,n)

f_fft = fft(f, N, 2);  % 2 indicates operation along columns

magnitude = abs(f_fft);

custom_map = [1 1 1; 1 0 0];  % [R G B] for white and red

figure (1);
imagesc(0:M-1, 0:N-1, magnitude);
colormap(custom_map); 
clim([0 128]);

xlabel('Sequence Index (m)');
ylabel('Frequency Bin (n)');
title('Spectrogram of 128-point FFT along Columns');

axis xy;