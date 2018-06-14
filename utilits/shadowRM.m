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
% 1. im_nf/flash: no_flash/flash image pair M*N  by  3
% 2. mask: mask for the image M by N
% 3. s_img: the image would like to remove the effect of flash shadows
% optional: 
% 4. imT: detected flash shadows
% 5. regions: the regions in the shadows are interested

%   

% Output
% 1. W1: the interpolated s_img for the detected shadow regions


function W1 = shadowRM(im_nf, im_f, mask, s_img, imT, regions)
    K = 50;
    [M, N, C] = size(im_nf);
    if nargin < 5
        imT = detectFlashShadows(im_nf, im_f, mask);
        %% if certain regions are given by the users.

        imT = reshape(imT, [M, N]);
        temp = zeros(size(imT));
        temp(regions.y, regions.x) = imT(regions.y, regions.x);
        imT = temp;

    end

    imT = imT(:);
    %%
    im_nf = reshape(im_nf, [], C);
    im_f = reshape(im_f, [], C);

    s_img = reshape(s_img, [], C);

    W1 = s_img;
    idS = find(imT == 1);
    [nf_mag,~] = imgradient(rgb2gray(reshape(im_nf, [M, N, C])));
    [f_mag,~] = imgradient(rgb2gray(reshape(im_f, [M, N, C])));

    %% filter out penumbra pixels
    % shadowMask = ones(M*N, 1) - imT;
    wSize = 3;%round(.01 * sqrt(M^2 + N^2)) + 1;
    hWinSize = 1;%(wSize - 1)/2;
    [row_idx, col_idx] = ind2sub([M, N], idS);

    penIDmat = zeros(M, N);
    for kk = 1:length(idS)
        [tempR,tempC] = meshgrid(row_idx(kk)-hWinSize:row_idx(kk)+hWinSize, ...
                        col_idx(kk)-hWinSize:col_idx(kk)+hWinSize); 
        tempR = tempR(:);
        tempC = tempC(:);
        tempR(floor(wSize/2) + 1) =[];
        tempC(floor(wSize/2) + 1) =[];

        % discard the out of boundary pixels
        discard_id = find(tempR <= 0 | tempC <= 0 | tempR > M | tempC > N); 

        tempR(discard_id) = [];
        tempC(discard_id) = [];
        % discard pixels out of the mask
        tempID = sub2ind([M N], tempR, tempC);
        discard_id = find(imT(tempID) == 1);
        tempID(discard_id) = [];

        penIDmat(tempID) = penIDmat(tempID) + 1;
    end
    imG = zeros(M * N, 1);
    imG(f_mag - nf_mag > 0.01) = 1;
    idTotal = find( ((penIDmat(:).*imG) > 0) & (imT == 0)); 
    idN = unique([idTotal;idS]);
    %imR = zeros(M, N, 1); imR(idN) = 1; %imshow(reshape(imR, [M, N]))
    fprintf('Begin to process %d shadow pixels.....\n', length(idN));

    
    s_idd = idN;% BW1 == 1; %

    tmp = zeros(M*N, 1);
    tmp(s_idd) = 1;
    us_idd = find(tmp == 0);


    [srow_id, scol_id] = ind2sub([M, N], idTotal);
    [usrow_id, uscol_id] = ind2sub([M, N], us_idd);

    
    sW = im_nf(idTotal, :);
    usW = im_nf(us_idd, :);
    repT  = length(us_idd);
    % id = [];
    fprintf('Process penumbra pixels.....\n');
    for kk = 1:length(srow_id)
       
        curr = repmat(sW(kk, :), repT, 1);
        dif_int = usW - curr;
        sumD = sqrt(dif_int(:, 1).^2  + dif_int(:, 2).^2 + dif_int(:, 3).^2); %
        sumD = sumD./sqrt(sum(sW(kk, :).^2));
        sumD = (sumD - min(sumD))./(max(sumD) - min(sumD));
        dif_row = usrow_id - srow_id(kk);
        dif_col = uscol_id - scol_id(kk);
        sumT = sqrt(dif_row.^2 + dif_col.^2);
        sumT = (sumT - min(sumT))./(max(sumT) - min(sumT));
        sumD = sumD+sumT;
        [~, id] = sort(sumD, 'ascend');
        ww = 1 - sumD(id(1:K))./sum(sumD(id(1:K)));
        ww = ww./sum(ww);
        W1(idTotal(kk), :) =  sum(s_img(us_idd(id(1:K)), :).*repmat(ww, [1 3]), 1);      
    end
    fprintf('Done.....\n');
    fprintf('Process umbra pixels.....\n');
    tmp = zeros(M*N, 1);
    tmp(idS) = 1;
    us_idd = find(tmp == 0);
    [srow_id, scol_id] = ind2sub([M, N], idS);
    [usrow_id, uscol_id] = ind2sub([M, N], us_idd);
    sW = im_nf(idS, :);
    usW = im_nf(us_idd, :);
    repT  = length(us_idd);
    % id = [];
    for kk = 1:length(srow_id)
        curr = repmat(sW(kk, :), repT, 1);
        dif_int = usW - curr;
        sumD = sqrt(dif_int(:, 1).^2  + dif_int(:, 2).^2 + dif_int(:, 3).^2); %
        sumD = sumD./sqrt(sum(sW(kk, :).^2));
        sumD = (sumD - min(sumD))./(max(sumD) - min(sumD));
        dif_row = usrow_id - srow_id(kk);
        dif_col = uscol_id - scol_id(kk);
        sumT = sqrt(dif_row.^2 + dif_col.^2);
        sumT = (sumT - min(sumT))./(max(sumT) - min(sumT));
        sumD = sumD+sumT;
        [~, id] = sort(sumD, 'ascend');
        ww = 1 - sumD(id(1:K))./sum(sumD(id(1:K)));
        ww = ww./sum(ww);
        W1(idS(kk), :) =  sum(s_img(us_idd(id(1:K)), :).*repmat(ww, [1 3]), 1);

    end
    fprintf('Done.....\n');
%     W = W1;%[WR WG WB];
%     W(isnan(W(:, 1)), 1) = 1;
%     W(isnan(W(:, 2)), 2) = 1;
%     W(isnan(W(:, 3)), 3) = 1;
%     W2 = reshape(W, [M, N, C]);
%     smoothness = 0.0001;%0.0001;
%     im_nf = reshape(im_nf, [M, N, C]);
%     for i = 1:3
%         W2(:,:,i) = imguidedfilter(W2(:,:,i), im_nf, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
%     end
%     W2 = W2./ repmat(sum(W2, 3), [1 1 3]);


end