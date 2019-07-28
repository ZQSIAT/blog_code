% by zqs
% turn picture
% 2018-3-22


% interpolation seting,if set '0' it is nearest neighbor interpolation,if set '1' it is bilinear interpolation,
% if set '2' it is bicubic interpolation,if set '3' it is area-based (or super) interpolation,
% if set '4' it is Lanczos interpolation over 8x8 neighborhood.

clc;
clear;
close all;
num_bins = 5;
num_methods = 3;

I_0 = load('Chapter2_1_interpolation_0.txt');
I_1 = load('turn_picture_interpolation_0.txt');
r_0 = corrcoef(I_0,I_1);
[psnr_0,mse_0] = psnr_mse(I_0,I_1);

I_2 = load('Chapter2_1_interpolation_1.txt');
I_3 = load('turn_picture_interpolation_1.txt');
r_1 = corrcoef(I_2,I_3);
[psnr_1,mse_1] = psnr_mse(I_2,I_3);

I_4 = load('Chapter2_1_interpolation_2.txt');
I_5 = load('turn_picture_interpolation_2.txt');
r_2 = corrcoef(I_4,I_5);
[psnr_2,mse_2] = psnr_mse(I_4,I_5);

I_6 = load('Chapter2_1_interpolation_3.txt');
I_7 = load('turn_picture_interpolation_3.txt');
r_3 = corrcoef(I_6,I_7);
[psnr_3,mse_3] = psnr_mse(I_6,I_7);

I_8 = load('Chapter2_1_interpolation_4.txt');
I_9 = load('turn_picture_interpolation_4.txt');
r_4 = corrcoef(I_8,I_9);
[psnr_4,mse_4] = psnr_mse(I_8,I_9);



% bins data

data = zeros(num_bins, num_methods);

data(1,1) = psnr_0;
data(2,1) = psnr_1;
data(3,1) = psnr_2;
data(4,1) = psnr_3;
data(5,1) = psnr_4;

data(1,2) = mse_0;
data(2,2) = mse_1;
data(3,2) = mse_2;
data(4,2) = mse_3;
data(5,2) = mse_4;

data(1,3) = r_0(1,2)*100;
data(2,3) = r_1(1,2)*100;
data(3,3) = r_2(1,2)*100;
data(4,3) = r_3(1,2)*100;
data(5,3) = r_4(1,2)*100;

% draw the graph
figure;
handle_bar = bar(data);
ax = gca;

title('Turn Picture Result');
ax.YGrid = 'on';
ax.XGrid = 'on';

ylabel('value');
xlabel('group');

leg_handle = legend('\bf{}PSNR', '\bf{}MSE', '\bf{}CORRELATION', 'Location', 'northeast');

% figure position and size : [left bottom width height]
set(gcf, 'Position', [100, 500, 400, 260]);

name = {'nearest', 'bilinear', 'bicubic', 'area-based', 'Lanczos'};

set(gca, 'XTickLabel', name);
% Place text atop the bar


function [PSNR,MSE] = psnr_mse(X,Y)
if any(size(X)~=size(Y))
   error('The input size is not equal to each other!');
end
D = X-Y;
MSE = sum(D(:).*D(:))/numel(X);
PSNR = 10*log10(255^2/MSE);
end