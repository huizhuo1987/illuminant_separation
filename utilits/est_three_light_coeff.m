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
% 1. gamma: M*N  by  3
% 2. mask: mask for the image M by N
% 3. cfactor: the threshold to find the corner points

%   

% Output
% 1. illum1, illum2, illum3: estimated lighting coefficents

function [illum1, illum2, illum3] = est_three_light_coeff(gamma, mask, cfactor)

    gamma = gamma(:, mask > 0);
    idM = mask > 0;

    center_gamma = mean(gamma, 2);
    center_gamma = center_gamma/norm(center_gamma);

    R = rot_world2local(center_gamma);
    gamma_proj = R*gamma;
    
%     alpha = 1./gamma_proj(3, :);
%     gamma_proj = gamma_proj.*repmat(alpha, [3 1]);

    qq = gamma_proj(1:2, :);
    [qq_density, density]= calDensityMap(qq, 200);
    
    idd = find(qq_density > cfactor);
    xScatter = qq(1, idd);
    yScatter = qq(2, idd);
    [trix, triy] = minboundtri(xScatter, yScatter, 0.1); %(xScatter, yScatter);
    g1 = [trix(1); triy(1)];
    g2 = [trix(2); triy(2)];
    g3 = [trix(3); triy(3)];
    

    [~, idd] = min( sum((gamma_proj(1:2, :) - repmat(g1(:), [1, size(gamma_proj, 2)])).^2, 1));
    g1 = gamma_proj(:, idd);
    [~, idd] = min( sum((gamma_proj(1:2, :) - repmat(g2(:), [1, size(gamma_proj, 2)])).^2, 1));
    g2 = gamma_proj(:, idd);
    [~, idd] = min( sum((gamma_proj(1:2, :) - repmat(g3(:), [1, size(gamma_proj, 2)])).^2, 1));
    g3 = gamma_proj(:, idd);
    

    illum1 = R'*g1;
    illum2 = R'*g2;
    illum3 = R'*g3;