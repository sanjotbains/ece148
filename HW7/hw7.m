set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');

% HW7.m - Estimation of Range Profiles

%% Range Estimation
%{
    The purpose of this exercise is to perform step-frequency FMCW radar 
    imaging using FFT for range estimation. 
 
    The data set was taken over a section of the walkway pavement in front 
    of the Broida Hall. The ground-penetrating radar imaging unit scanned 
    along a linear path and took data at 200 spatial positions. The spatial 
    spacing between the data-collection positions is 0.0213 m (2.13 cm). 
 
    At each data-collection position, the system illuminates the subsurface 
    region of the walkway area with microwaves in the step-frequency (FMCW) 
    mode, stepping through 128 frequencies with a constant increment, from 
    0.976 GHz to 2.00 GHz. The relative permittivity εr is approximately 6.0.  
 
    The data set is in the form of a (200 x 128) array, corresponding to 128
    complex-amplitude values (for the 128 frequencies) at 200 data-
    acquisition positions. 
 
    Your task is to perform the image reconstruction of the subsurface 
    profile. The result will be an image compiled in the form of 200 depth 
    profiles corresponding to the 200 data-acquisition positions.
%}

% Load the data
gpr_data = load('gpr_data/gpr_data.mat');

F = gpr_data.F; 
    % The 128x200 data array of FMCW frequency data collected during the
	% experiment (where each column corresponds to the antenna position,
	% and each row corresponds to a frequency value). 
f = gpr_data.f;
    % 128 point frequency vector corresponding to the transmit frequency 
    % band used in (Hz).
da = gpr_data.da;
    % Position of the antenna at each data collection point (m).


%% RANGE ESTIMATION

% (a) The propogation speed needs to be adjusted by the relative permittivity.

epsilon_r = 6;

%{
    The propagation speed of electromagnetic waves in a dielectric medium is 
    inversely proportional to the square root of the relative permittivity. 
    The formula for the propagation speed (v) is v = c / sqrt(epsilon_r), 
    where c is the speed of light in a vacuum.
%}

c = 3 * 10^8; # m/s 
velocity = c / sqrt(epsilon_r); # = 1.2247e+08 m/s

bandwidth = 1 * 10^9; % 1 GHz

% Time Resolution
delta_t = 1 / bandwidth; % 10^-9 s

% (b) The depth profile needs to be scaled by a factor of two for the round-
%     trip progation.

% Range Resolution
delta_r = velocity * delta_t / 2; % = 0.061237 m

% Frequency Resolution
delta_f = (2 - 0.976) * 10^9 / 128; % = 8 MHz

% Max Time 
range_time_max = 1 / delta_f; % = 1.25 * 10^-7 s

% Max Distance
range_distance_max = velocity * range_time_max / 2; % = 7.6547 m

% (c) The intensity of the depth profiles near the transceiver is due to strong 
%     surface reflection. It should be gated out in order to visualize the 
%     depth profiles adequately. 


% (d) The rebars are at shallow depths. Thus, the displayed profiles up to 
%     20-30 cm will be sufficient.

% Dividing 765 cm / 24 cm (within our range of interest), yields approximately
% 32.

% We have 128 frequencies: 128/32 = 4, meaning we only care about the first 4 px.

% Show the original data, only the first 8 pixels
figure (1);
imagesc( abs(
    fft(F)(121:128, :) 
    ));    
colormap('turbo');
title('Range Profiles');


%% INTERPOLATION

%% Interpolating Range (Depth) for Display (Columns)

% THE SIZE OF OUR DATA IS 128 columns x 200 rows.

% We want to interpolate the depth data to 2048 points.
% Interpolation factor = 2048 / 128 = 16

% We need to zero-pad the data to 2048 points.
% Very naively, we can just pad with zeros at the end of each column.
% We know the freq resolution is 8 Mhz, and the minimum frequency is 0.976 GHz.
% That means we can determine the number of empty points we should add infront
% of the data.

% Create the padded data
% We are padding the frequency axis
pts_front = 0.976 * 10^9 / delta_f; % = 122
F_padded = zeros(2048, 200);
F_padded(pts_front + 1 : pts_front + 128, :) = F;

% THE SIZE OF OUR DATA IS 2048 rows x 200 columns.

% Perform the IFFT on the padded data
% The IFFT is performed along the columns.

F_ifft = ifft(F_padded, [], 1);

% The IFFT turns the columnar frequency data into depth profiles.

F_ifft = ifft(F_padded, [], 1)(1985:2048, :);
    % Number of pixels for display = (4*16) = 64
    % Pixel size = 6.1237 cm / 16 = 0.3827 cm

% THE SIZE OF OUR DATA IS 64 rows x 200 columns.

figure (2);
imagesc(abs(F_ifft));
colormap('turbo');
title('Range Profiles');
xlabel('Antenna Position (m)');
ylabel('Depth (m)');
print('range_interp.png', '-dpng', '-tight', '-r300');

%% Interpolating receiver positions. (horizontal)

% We have 200 receiver positions

% We will first pad by 28 zeros at the beginning and end of the receiver pos axis.
% This lets us FFT the data.
F_padded = zeros(64, 256);
F_padded(:, 29:228) = F_ifft;

# THE SIZE OF OUR DATA IS 64 rows x 256 columns.

% Now we can perform the FFT along the position axis.
F_fft = fft(F_padded, [], 2);

% In the frequency domain we can zero pad the position data to 2048 points.
% We add zeros to the middle of the data
F_padded = zeros(64, 2048);
F_padded(:, 1:128) = F_fft(:, 1:128);
F_padded(:, 1921:2048) = F_fft(:, 129:256);

# THE SIZE OF OUR DATA IS 64 rows x 2048 columns.

% Perform the IFFT on the padded data
F_ifft = ifft(F_padded, [], 2);

% We are now back in the spatial domain, with 2048 receiver points.
% The interpolation ratio is 2048/256 = 8.
% The pixel size is 0.0213 m / 8 = 0.0026625 m = 2.6625 mm = 0.26625 cm

% To get a properly scaled image we must take into account the disparate
% pixel sizes in the depth and horizontal directions.

% Each depth pixel is 0.3827 cm, and each horizontal pixel is 0.26625 cm.
% We can use the `imagesc` function to display the data, and set the
% 'XData' and 'YData' properties to scale the image correctly.

figure (3);
imagesc(abs(F_ifft)(:, 225:1825));
colormap('jet');
set(gca, 'YDir', 'normal'); % Set Y-axis direction to normal
% set(gca, 'XTick', 1:8:2048, 'XTickLabel', round((da(1:8:end) * 100) / 100));
% xlabel('Antenna Position (m)');
% ylabel('Depth (m)');
% title('Range Profiles with Interpolated Receiver Positions');
axis equal; % Set equal aspect ratio
% set(gca, 'DataAspectRatio', [0.26625 0.3827 1]); % Set aspect ratio
print('range_interp_horizontal.png', '-dpng', '-tight', '-r300');

img = abs(F_ifft)(:, 225:1825);
img = img - min(img(:));
img = img / max(img(:));
img_uint8 = uint8(img * 255);
cmap = jet(256);
imwrite(img_uint8, cmap, 'depth_map.png');


