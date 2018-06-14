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
% 1. L: measured illuminant spectra data N by M, 
%    N: the samples for the spectra, M the number of samples in the database   
% 2. R: measured reflectance spectra data N by M
% 3. C: camera response spectra
% 4. option: the difference methods in generating the basis

%   

% Output
% 1. UR, UL: the generated basis for reflectance and illumination,
%    respectively



function [UR, UL] = genBase(L, C, R, option)

switch option
    case 'joint'
        [VR, ~, VL] = svds(R'*diag(C*ones(3,1))*L, 3);
        UR = orth(R*VR);
        UL = orth(L*VL);
          
    case 'pca'

        [UR, SigmaR, VR] = svds(R, 3);
        
        [UL, SigmaL, VL] = svds(L, 3);
        UR = orth(UR);
        UL = orth(UL);
    case 'wpca'

        R_total = [R.*repmat((C(:, 1)), 1, size(R, 2)) R.*repmat((C(:, 2)), 1, size(R, 2)) R.*repmat((C(:, 3)), 1, size(R, 2))];
        L_total = [L.*repmat((C(:, 1)), 1, size(L, 2)) L.*repmat((C(:, 2)), 1, size(L, 2)) L.*repmat((C(:, 3)), 1, size(L, 2))];
        [UR, ~, ~] = svds(R_total, 3);       
        [UL, ~, ~] = svds(L_total, 3);   

        
        UR = orth(UR);
        UL = orth(UL);
    case 'wpca1'
        
        [UR(:, 1),~,~] = svds(R.*repmat((C(:, 1)), 1, size(R, 2)), 1);
        [UR(:, 2),~,~] = svds(R.*repmat((C(:, 2)), 1, size(R, 2)), 1);
        [UR(:, 3),~,~] = svds(R.*repmat((C(:, 3)), 1, size(R, 2)), 1);

        [UL(:, 1),~,~] = svds(L.*repmat((C(:, 1)), 1, size(L, 2)), 1);
        [UL(:, 2),~,~] = svds(L.*repmat((C(:, 2)), 1, size(L, 2)), 1);
        [UL(:, 3),~,~] = svds(L.*repmat((C(:, 3)), 1, size(L, 2)), 1);
        
        UR = orth(UR);
        UL = orth(UL);
end


end