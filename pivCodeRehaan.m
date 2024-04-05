clear all
close all
clc

%Begin with first two frames

%Read image pair 
image_a = im2double(imread('500fps_1000ss_rehaan_maggie_C001H001S0001000050.bmp'));
image_b = im2double(imread('500fps_1000ss_rehaan_maggie_C001H001S0001000051.bmp'));

image_dimensions = size(image_a); %in pixels
image_size = image_dimensions(1); %assuming square image

%Subtract average intensity
image_a = image_a - mean(image_a,'all'); 
image_b = image_b - mean(image_b,'all');

%Select IA size (recommended 64x64 or 48x48)
IA_size = 64; %in pixels
IA_overlap = IA_size*0.25; 

%Setting floor to zero
image_a = image_a - min(image_a,[],'all'); 
image_b = image_b - min(image_b,[],'all');

%Using a Hamming window to eliminate end effects
window = hamming(IA_size)*hamming(IA_size)';

%{
For showing preprocessing difference in report (not applying Hamming 
window)...uncomment to run
%}
% window = IA_size*IA_size;

%Initializing matrices
x_vel = zeros(size(28,61));
y_vel = zeros(size(28,61));
correlation_value = zeros(size(28,61));
velocity_magnitude = zeros(size(correlation_value));
correlation_ratios = zeros(size(27,61));

%Looping over the IAs find cross correlation between the image pairs
for i = 1:(image_size-IA_size)/IA_overlap+1
    for j = 1:(image_size-IA_size)/IA_overlap+1
       
        window_a = image_a(1+(i-1)*IA_overlap: (i-1)*IA_overlap+IA_size, 1+(j-1)*IA_overlap: (j-1)*IA_overlap+IA_size).*window;
        window_b = image_b(1+(i-1)*IA_overlap: (i-1)*IA_overlap+IA_size, 1+(j-1)*IA_overlap: (j-1)*IA_overlap+IA_size).*window;
        correlation = xcorr2(window_a, window_b);
        
        %Normalizing the correlation
        %First find autocorr of each IA
        autocorr_a = xcorr2(window_a,window_a);
        autocorr_a = autocorr_a(IA_size,IA_size); %use center value
        autocorr_b = xcorr2(window_b,window_b);
        autocorr_b = autocorr_b(IA_size,IA_size);
        %Normalize using sqrt of the autocorrelations
        normalize = sqrt(autocorr_a*autocorr_b);
        norm_correlation = correlation/normalize;
        [x, y] = ndgrid((-IA_size + 1):(IA_size - 1));
        interpolation = griddedInterpolant(x, y, norm_correlation,'spline');
        [peaks, value] = fminsearch(@(f) -interpolation(f(1), f(2)), [0, 0]);

        x_vel(i, j) = peaks(2);
        y_vel(i, j) = peaks(1);
        correlation_value(i, j) = -value; 
        velocity_magnitude(i,j) = sqrt(peaks(2)^2 + peaks(1)^2);

        %Function below found from this forum
        % https://uk.mathworks.com/matlabcentral/answers/277512-how-to-find-peaks-in-3d-mesh
       % ix = find(imregionalmax(norm_correlation));
       % list = sort(norm_correlation(ix), 'descend');

       % correlation_ratios(i,j) = list(1)/list(2);

    end
end


for i = 1:(image_size-IA_size)/IA_overlap+1
    for j = 1:(image_size-IA_size)/IA_overlap+1
        x_pos(i, j) = j;
        y_pos(j, i) = j;
    end
end

figure
quiver(x_pos, y_pos, x_vel, y_vel)
axis ij
title('PIV Relative Velocity Vector Field')
xlabel('x-position')
ylabel('y-position')

figure 
surf(x_pos, y_pos, curl(x_pos, y_pos, x_vel, y_vel))
view(2)
colorbar
title('PIV Relative Vorticity Map')
xlabel('x-position')
ylabel('y-position')
xlim([0,32])
ylim([0,32])

axis ij

figure 
surf(x_pos, y_pos, correlation_value)
view(2)
colorbar
title('Max Correlation Map')
xlabel('x-position')
ylabel('y-position')
axis ij
xlim([0,32])
ylim([0,32])

figure 
surf(x_pos, y_pos, velocity_magnitude)
view(2)
colorbar
title('Velicty Magnitude Plot')
xlabel('x-position')
ylabel('y-position')
axis ij
xlim([0,32])
ylim([0,32])

%% Run all image pairs

parameter = "Enter '1' to run all image pairs: ";
run = input(parameter);

if run==1
folder = '/Users/mageean/Documents/MATLAB/experimental-fluids-PIV/experimental_data';
f =dir([folder '/*.bmp']);
n = length(f); %calculate number of files
pairs = n/2;

%Loop over all image pairs
for i =1:pairs
end
%Read image pair
image_a = im2double(imread('c000a.bmp'));
image_b = im2double(imread('c000b.bmp'));

%Subtract average intensity
image_a = image_a - mean(image_a(:),'all'); 
image_b = image_b - mean(image_b(:),'all');

%Select IA size (recommended 64x64 or 48x48)
image_size = sqrt(numel(image_a)); %in pixels, assuming square image
IA_size = 48; %in pixels
IA_overlap = IA_size*0.25; 

%Using a Hamming window to eliminate end effects
window = hamming(IA_size)*hamming(IA_size)';

for i = 1:(image_size-IA_size)/IA_overlap+1
    for j = 1:(image_size-IA_size)/IA_overlap+1
       
        window_a = image_a(1+(i-1)*IA_overlap: (i-1)*IA_overlap+IA_size, 1+(j-1)*IA_overlap: (j-1)*IA_overlap+IA_size).*window;
        window_b = image_b(1+(i-1)*IA_overlap: (i-1)*IA_overlap+IA_size, 1+(j-1)*IA_overlap: (j-1)*IA_overlap+IA_size).*window;
        norm_correlation = normxcorr2(window_a, window_b);
  
        [x, y] = ndgrid((-IA_size + 1):(IA_size - 1));
        interpolation = griddedInterpolant(x, y, norm_correlation,'spline');
        [peaks, value] = fminsearch(@(f) -interpolation(f(1), f(2)), [0, 0]);

        x_vel(i, j) = peaks(2);
        y_vel(i, j) = peaks(1);
        correlation_value(i, j) = -value; 
        velocity_magnitude(i,j) = sqrt(peaks(2)^2 + peaks(1)^2);

        %Function below found from this forum
        % https://uk.mathworks.com/matlabcentral/answers/277512-how-to-find-peaks-in-3d-mesh
        ix = find(imregionalmax(norm_correlation));
        list = sort(norm_correlation(ix), 'descend');

        correlation_ratios(i,j) = list(1)/list(2);

    end
end


for i = 1:(image_size-IA_size)/IA_overlap+1
    for j = 1:(image_size-IA_size)/IA_overlap+1
        x_pos(i, j) = j;
        y_pos(j, i) = j;
    end
end

figure
quiver(x_pos, y_pos, x_vel, y_vel)
axis ij
title('PIV Relative Velocity Vector Field')
xlabel('x-position')
ylabel('y-position')

figure 
surf(x_pos, y_pos, curl(x_pos, y_pos, x_vel, y_vel))
view(2)
colorbar
title('PIV Relative Vorticity Map')
xlabel('x-position')
ylabel('y-position')
xlim([0,32])
ylim([0,32])

axis ij

figure 
surf(x_pos, y_pos, correlation_value)
view(2)
colorbar
title('Max Correlation Map')
xlabel('x-position')
ylabel('y-position')
axis ij
xlim([0,32])
ylim([0,32])

figure 
surf(x_pos, y_pos, velocity_magnitude)
view(2)
colorbar
title('Velicty Magnitude Plot')
xlabel('x-position')
ylabel('y-position')
axis ij
xlim([0,32])
ylim([0,32])
else 
end