%This code breaks down images into N x N interrogation areas (IAs) and cross-
%correlates these areas between two images in order to calculate velocity
%and vorticity fields.
%This will also include windowing and overlapping IAs.
%Next, we will estimate the quality of the calculated velocity vector
%written by Raymond Gresalfi, Adam Poche, and Maggie Brown
%ENGN2912T: Experimental Fluid Mechanics - HW3
%Spring 2024 - Professor Breuer, Brown University

%MAKE SURE YOU HAVE SIGNAL AND IMAGE PROCESSING TOOLBOXES TO RUN
% -----------------------------------------------------------------------
%--------------------------- SET UP ---------------------------
clear all
clc;
%Reach each frame, remove mean and min
im1raw = imread("zoe_maggie_piv_500fps_2000ss_04_02_2024_1_C001H001S0001000303.bmp");
im1sub = im1raw - mean(im1raw,'all')-min(min(im1raw));
im2raw = imread("zoe_maggie_piv_500fps_2000ss_04_02_2024_1_C001H001S0001000305.bmp");
im2sub = im2raw - mean(im2raw,'all')-min(min(im2raw));
%Hanning Window
[m,n,~] = size(im1sub);
hw = hanning(n) * hanning(m)';
im1 = im2double(im1sub) .* hw;
im2 = im2double(im2sub) .* hw;
%Set Scales
delta_t = 1; % seconds
pix_scale = 1; %meters/pixel
edge_pixel = length(im1); %assume size of each frame is constant/square
%Establish Interrogation Areas + Overlap
N = 16; %size of interrogation areas
overlap = 4; %number of overlaps per IA
I = edge_pixel/N; %number of interrogation areas (before overlap)
IA_num = (I-1)*overlap + 1; %number of interrogation areas (after overlap)
%Create stored vectors for loop
x = zeros(1,IA_num-1);
y = zeros(1,IA_num-1);
u = zeros(IA_num-1);
v = zeros(IA_num-1);
vel = zeros(IA_num-1);
%Quality vectors
ratio = zeros(IA_num-1);
max_corr = zeros(IA_num-1);
%--------------------------- LOOP ---------------------------
%Loop over all interrogation areas
for i = 1:IA_num-1
    for j = 1:IA_num-1
        %Reset maximum
maximum = 0;
        %Define bounds for each overlapped area
        leftbound = (N/overlap)*(i-1)+1;
        rightbound = leftbound + N;
        topbound = (N/overlap)*(j-1)+1;
        botbound = topbound + N;
        IA1 = im1(leftbound:rightbound,topbound:botbound);
        IA2 = im2(leftbound:rightbound,topbound:botbound);
        %Skip iteration if image is all black
if IA1 == 0
    continue
end

if IA2 == 0
    continue 
end
        %Cross Correlate - Find Location of the Maximums
        cross_corr = normxcorr2(IA2,IA1);
        maximum = max(max(cross_corr));
        [un,in] = unique(cross_corr);
        second_max = un(end-1);
        third_max = un(end-2);
        max_corr(i,j) = maximum;
        ratio(i,j) = maximum/second_max;
        %Skip iteration if peak does not look good!
        if maximum < 0.50
            continue
        end
        if abs(maximum - second_max) < 0.10
continue 
        end
if abs(third_max) > 0.70
   
   continue 
end
        %Find the maximum location
        [y_corr,x_corr] = find(cross_corr == maximum);
        %Skip iteration if two exact maximums
        if length(y_corr) ~= 1
            continue
        end
        if length(x_corr) ~= 1
            continue 
        end
       %Skip iteration if maximum too close to edge
        if abs(x_corr - 2*N) < 2
            continue
        end
        if abs(y_corr - 2*N) < 2
continue 
        end
        if x_corr <= 2
            continue
        end
        if y_corr <=2
continue 
        end
        %Manual 2D Gaussian Fit within small window
        reach = 2;
        resolution = 100;
        y_strip = cross_corr(y_corr,x_corr-reach:x_corr+reach);
        x_strip = cross_corr(y_corr-reach:y_corr+reach,x_corr);
        y_strip_coords = 1:1:length(y_strip);
        x_strip_coords = 1:1:length(x_strip);
        y_gauss_coords = 1:(1/resolution):length(y_strip);
        x_gauss_coords = 1:(1/resolution):length(x_strip);
        
        y_gauss = fit(y_strip_coords.',y_strip.','gauss1');
        x_gauss = fit(x_strip_coords.',x_strip,'gauss1');

        [max_gauss_y,iy] = max(y_gauss(y_gauss_coords));
        [max_gauss_x,ix] = max(x_gauss(x_gauss_coords));
        %Move correlation axes to middle
        x_corr_norm = (x_corr - size(IA1,1)) + (iy-(reach*resolution))/resolution; %gaussian correction
        y_corr_norm = y_corr - size(IA1,1) + (ix-(reach*resolution))/resolution;
        %Skip iteration if velocity oddly high
   
if abs(sqrt(x_corr_norm^2+y_corr_norm^2)) > 2
  continue
end
        %Calculate Velocities with Scale
        u(i,j) = x_corr_norm*pix_scale/delta_t;
        v(i,j) = y_corr_norm*pix_scale/delta_t;
        vel(i,j) = sqrt(u(i,j)^2+v(i,j)^2);
        %Establish quiver grid
        x(i) = i*N/overlap;
        y(j) = j*N/overlap;
        %~~~~~Insider test!~~~~~~
        if i == 44 && j == 44
            cross_corr_test = cross_corr;
            x_fit = 1:1:length(cross_corr_test);
            y_fit = x_fit;
            cross_corr_test_small = cross_corr_small;
            figure(1)
            surf(cross_corr)
            title('Test Correlation')
            disp('--- Random Results ---');
            disp(['maximum correlation = ',num2str(maximum)]);
            disp(['x disp = ',num2str(x_corr_norm),' pixels']);
            disp(['y disp = ',num2str(y_corr_norm),' pixels'])
            disp(['u = ',num2str(u(i,j)),' pixels/s']);
            disp(['v = ',num2str(v(i,j)),' pixels/s']);
end
        %~~~~~~~~~~~~~~~~~~~~~~~~
    end 
end
circulation = zeros(length(x));
vorticity = zeros(length(x));
for i = 2:length(x)-1
    for j = 2:length(y)-1
        %Calculate Vorticity using Circulation Integral
        circulation(i,j) = ...
           -u(i-1,j+1)- u(i,j+1)    -u(i+1,j+1) ...
           -v(i-1,j+1)              +v(i+1,j-1) ...
           -v(i-1,j)                +v(i+1,j)- ...
           -v(i-1,j-1)              +v(i+1,j-1) ...
           +u(i-1,j-1) +u(i,j-1)    +u(i+1,j-1);
        
        vorticity(i,j) = circulation(i,j)/9;
    end 
end


   %--------------------------- PLOT RESULTS ---------------------------
%PIV Overlay
figure(2)
imshow(im1)
hold on
figure(2)
quiver(x,y,u,-v,3,'g') %reverse v to fit image plot
hold off
%Velocity
figure(3)
imagesc(vel)
colorbar
title('Velocity')
%Vorticity
figure(4)
imagesc(vorticity)
colorbar
title('Vorticity')
%Highest Correlation Peak
figure(5)
imagesc(max_corr)
colorbar
title('Highest Correlation Peak')
%Ratio R1/R2
figure(6)
imagesc(1./ratio)
colorbar
title('Ratio R2/R1')