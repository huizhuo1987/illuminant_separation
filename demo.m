% Code for Illuminant Spectra-based Source Separation Using Flash Photography
% This code is based on the algorithm proposed in the paper
% "Illuminant Spectra-based Source Separation Using Flash Photographye", CVPR 2018
% Zhuo Hui, Kalyan Sunkavalli, Sunil Hadap, Aswin C. Sankaranarayanan
% When you use the code to build your algorithm, please cite this paper. 
% 
% Please contact the author Zhuo Hui if you have any problems with the code
% huizhuo1987@gmail.com
% 
% Copy rights reserved by the authors listed above.

% set up the route
setup()

% load the data to process
disp('Loading data...')
load('data\reflectance_illum_camera.mat')
load('data\images\two_lights\book.mat')
disp('Done.')

% generate basis for both reflectance/illumination
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Generate basis....')
[UR, UL] = genBase(L, C, R, 'wpca1');
for ind = 1:3
    E(:,:,ind) = UR'*diag(C(:, ind))*UL;
end
opt.E = E;

flash_light = .025*ones(size(L(:, 1)));
f = UL'*flash_light;
f = f/norm(f);
opt.f = f;

%% Demo for TWO lights

%% These parameter can be adjusted for better performance
opt.lambda = 1e-5;
opt.cutoff = .1;
opt.color_correct = 1;
opt.light_number = 2;
opt.shadow_mask = 1 - mask;
disp('Done.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%
disp('Processing the data....')
results = solve_light_sep(im_nf, im_f, mask, opt);
disp('Done');

%% Visulize the results
imshow([results.im1 results.im2].^(1/2.2))


%% Demo for Three lights

%% These parameter can be adjusted for better performance
disp('Loading data...')
load('data\images\three_lights\toy.mat')
disp('Done.')

opt.lambda = 1e-5;
opt.cutoff = .1;
opt.color_correct = 1;
opt.light_number = 3;
opt.shadow_mask = shadow_mask;
disp('Done.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%
disp('Processing the data....')
results = solve_light_sep(im_nf, im_f, mask, opt);
disp('Done');

%% Visulize the results
imshow([results.im1 results.im2].^(1/2.2))