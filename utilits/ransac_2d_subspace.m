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

%% This function is to implement the RANSAC function


function [n0, nx] = ransac_2d_subspace(pts, inlier_angle, iter, inliers, M)
    
    if ~exist('M', 'var')
        M = 2;
    end
    if ~exist('inliers', 'var')
        inliers = floor(size(pts, 2)/2);
    end
    if ~exist('inlier_angle', 'var')
        inlier_angle = (1*pi/180);
    end
    if ~exist('iter', 'var')
        iter = 100;
    end

    max_inlier = 0;
    for kk=1:iter
        x0 = randperm(size(pts, 2), 2);
        n0 = cross(pts(:, x0(1)), pts(:, x0(2)) );
        n0 = n0/norm(n0);

        ang = abs(asin(n0'*pts));
        in_count = length(find(ang <= inlier_angle));
        if (max_inlier <= in_count)
            max_inlier = in_count;
            n0_max = n0;

        end

    end

    ang = abs(asin(n0_max'*pts));
    % ang = abs(acos(n0_max'*pts));
    in_idx = find(ang <= inlier_angle);

    fprintf('Found %1.3f % \n', max_inlier/size(pts, 2));

    [U, S, V] = svds(pts(:, in_idx), 3);
    n0 = U(:, 3);
    nx = U(:, 1:2);
