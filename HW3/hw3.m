img = imread('hw3.png'); % Read the image
img = rescale(im2double(img)); % Convert the image to double precision

%% Problem 1
%{
    (Impulse invariance) Formulate a 3 x 3 matrix operator as the moving average filter in 
    the discrete form and apply the moving average filter to the given image for smoothing. 
    Show the result of the lowpass-filtering process.
%}

movingAverageFilter = [1 1 1; 1 1 1; 1 1 1] / 9; % Normalized moving average filter

smoothedImg = rescale(conv2(img, movingAverageFilter, 'same')); % Apply the filter to the image

avgFilter = figure(1); % Create a new figure
imshow(smoothedImg); % Display the smoothed image
%title('Smoothed Image using Moving Average Filter'); % Title for the figure

exportgraphics(avgFilter, 'smoothed_image.png', 'Resolution', 300); % Save the figure as a PNG file

clear movingAverageFilter; % Clear the filter variable to free memory

%% Problem 2
%{
    Produce composite images to illustrate the smoothing effect in the form of weighted 
    superpositions onto the original image.
%}

weights = [0.0, 0.33, 0.67, 1.0]; % Define weights for blending

figure(2); % Create a new figure for composite images
smoothComp = tiledlayout(2,2); % Create a 2x2 tiled layout for subplots

for i = 1:length(weights)
    compositeImg = weights(i) * smoothedImg + (1 - weights(i)) * img; % Weighted superposition
    nexttile; % Move to the next tile in the layout
    imshow(compositeImg); % Display composite image
    title(['Weight = ', num2str(weights(i))]); % Title for each subplot
end

% title(smoothComp, 'Composite Images with Smoothing Effect'); % Title for the tiled layout

exportgraphics(smoothComp, 'smoothed_composites.png', 'Resolution', 300); % Save the figure as a PNG file

clear smoothedImg compositeImg; % Clear variables to free memory

%% Problem 3
%{
    (Equation conversion) Formulate a 3 x 3 matrix operator as the gradient filter in the 
    discrete form and apply the gradient filter to the image for edge detection. Your results 
    should include (a) edge profile of the images in the horizontal direction, (b) edge profile 
    in the vertical direction, and (c) the combined edge profile.
%}

gradientFilterX = [-1 0 1; -2 0 2; -1 0 1]; % Gradient filter for horizontal edges
gradientFilterY = [-1 -2 -1; 0 0 0; 1 2 1]; % Gradient filter for vertical edges
gradientFilterCombined = [(-1+1j) 2j (1+1j); -2 0 2; (-1-1j) -2j (1-1j)]; % Combined gradient filter

xDerivative = rescale(conv2(img, gradientFilterX, 'same')); % Apply horizontal gradient filter
yDerivative = rescale(conv2(img, gradientFilterY, 'same')); % Apply vertical gradient filter
combinedDerivative = rescale(abs(conv2(img, gradientFilterCombined, 'same'))); % Apply combined gradient filter

figure(3); % Create a new figure for edge profiles
gradFilters = tiledlayout(1,3); % Create a 3x1 tiled layout for subplots

nexttile; % Move to the next tile in the layout
imshow(xDerivative, []); % Display horizontal edge profile
title('Alex is Dumb Edge Profile'); % Title for the subplot

nexttile; % Move to the next tile in the layout
imshow(yDerivative, []); % Display vertical edge profile
title('Vertical Edge Profile'); % Title for the subplot

nexttile; % Move to the next tile in the layout
imshow(combinedDerivative, []); % Display combined edge profile
title('Combined Edge Profile'); % Title for the subplot

% title(gradFilters, 'Edge Profiles'); % Title for the tiled layout

exportgraphics(gradFilters, 'edge_profiles.png', 'Resolution', 300); % Save the figure as a PNG file

clear gradientFilterX gradientFilterY gradientFilterCombined; % Clear filter variables

%% Problem 4
%{
    Produce composite images to illustrate the edge-detection effect in the form of weighted 
    superpositions onto the original image, separated into 3 different 2x2 plots.
%}

% Horizontal Edge Composite Images
figure(4); % Create a new figure for horizontal edge composites
xComp = tiledlayout(2, 2); % Create a 2x2 tiled layout
for i = 1:length(weights)
    nexttile; % Move to the next tile in the layout
    compositeImgX = weights(i) * xDerivative + (1 - weights(i)) * img; % Weighted superposition for x-derivative
    imshow(compositeImgX, []); % Display composite image for x-derivative
    title(['Weight = ', num2str(weights(i)), ' (X)']); % Title for each subplot
end
% title(xComp, 'Composite Images with Horizontal Edge Detection'); % Title for the tiled layout

% Vertical Edge Composite Images
figure(5); % Create a new figure for vertical edge composites
yComp = tiledlayout(2, 2); % Create a 2x2 tiled layout
for i = 1:length(weights)
    nexttile; % Move to the next tile in the layout
    compositeImgY = weights(i) * yDerivative + (1 - weights(i)) * img; % Weighted superposition for y-derivative
    imshow(compositeImgY, []); % Display composite image for y-derivative
    title(['Weight = ', num2str(weights(i)), ' (Y)']); % Title for each subplot
end
% title(yComp, 'Composite Images with Vertical Edge Detection'); % Title for the tiled layout

% Combined Edge Composite Images
figure(6); % Create a new figure for combined edge composites
combinedComp = tiledlayout(2, 2); % Create a 2x2 tiled layout
for i = 1:length(weights)
    nexttile; % Move to the next tile in the layout
    compositeImgCombined = weights(i) * combinedDerivative + (1 - weights(i)) * img; % Weighted superposition for combined derivative
    imshow(compositeImgCombined, []); % Display composite image for combined derivative
    title(['Weight = ', num2str(weights(i)), ' (Combined)']); % Title for each subplot
end
% title(combinedComp, 'Composite Images with Combined Edge Detection'); % Title for the tiled layout

exportgraphics(xComp, 'horizontal_edge_composites.png', 'Resolution', 300); % Save the horizontal edge composites as a PNG file
exportgraphics(yComp, 'vertical_edge_composites.png', 'Resolution', 300); % Save the vertical edge composites as a PNG file
exportgraphics(combinedComp, 'combined_edge_composites.png', 'Resolution', 300); % Save the figure as a PNG file

clear xDerivative yDerivative combinedDerivative; % Clear variables to free memory
clear compositeImgX compositeImgY compositeImgCombined; % Clear composite image variables

%% Problem 5
%{ 
    (Equation conversion) Formulate a 3 x 3 matrix operator as the Laplacian filter in the 
    discrete form and apply the Laplacian filter to the given image for peak detection. Show 
    the distribution of peaks.
%}

laplacianFilter = [0 -1 0; -1 4 -1; 0 -1 0]; % Laplacian filter

laplacianImg = rescale(conv2(img, laplacianFilter, 'same')); % Apply Laplacian filter

laplacePlot = figure (7); % Create a new figure for peak detection
imshow(laplacianImg, []); % Display the Laplacian image
% title('Laplacian Filtered Image for Peak Detection'); % Title for the figure

exportgraphics(laplacePlot, 'laplacian_image.png', 'Resolution', 300); % Save the figure as a PNG file

clear laplacianFilter; % Clear the filter variable to free memory

%% Problem 6
%{
     Produce composite images to illustrate the peak-detection effect in the form of weighted 
    superpositions onto the original image. 
%}

figure(8); % Create a new figure for composite images
laplaceComp = tiledlayout(2,2); % Create a 2x2 tiled layout for subplots

for i = 1:length(weights)
    compositeImgLaplacian = weights(i) * laplacianImg + (1 - weights(i)) * img; % Weighted superposition for Laplacian
    nexttile; % Move to the next tile in the layout
    imshow(compositeImgLaplacian, []); % Display composite image for Laplacian
    title(['Weight = ', num2str(weights(i))]); % Title for each subplot
end

% title(laplaceComp, 'Composite Images with Peak Detection Effect'); % Title for the tiled layout

exportgraphics(laplaceComp, 'laplacian_composite_images.png', 'Resolution', 300); % Save the figure as a PNG file

clear laplacianImg compositeImgLaplacian; % Clear variables to free memory


