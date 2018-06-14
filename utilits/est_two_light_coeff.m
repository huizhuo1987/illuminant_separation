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

%% This function is to implement the esimtation of lighting coefficients

% Input
% 1. gamma: 3 by  M*N
% 2. mask_t: mask for the image M by N
% 3. cfactor: the threshold to find the corner points

%   

% Output
% 1. illum1, illum2: estimated lighting coefficents



function [illum1, illum2] = est_two_light_coeff(gamma, mask_t, cfactor)
    [n0, nx] = ransac_2d_subspace(gamma(:, mask_t>0), .3*pi/180, 1000);
    theta = (0:.1:360)*pi/180;
    circ = nx*[cos(theta); sin(theta)];

    iit = find((mask_t) > 0);
    for ii = 1:length(iit)
        [~, idx(ii)] = max(circ'*gamma(:, iit(ii)));
    end
    [h_idx] = hist(idx, 1:size(circ, 2));

    cutoff = cfactor*sum(h_idx)/length(theta);

    theta1 = theta(min(find(h_idx > cutoff)));
    theta2 = theta(max(find(h_idx > cutoff)));
    illum1 = nx*[cos(theta1); sin(theta1)];
    illum2 = nx*[cos(theta2); sin(theta2)];

