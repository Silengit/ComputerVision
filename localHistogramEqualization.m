function [img1] = localHistogramEqualization(img0)
%Local histogram equalization for 2D colorful image
rawSize = size(img0);
block = 4;
blockHeight = floor(rawSize(1)/block);
blockWidth = floor(rawSize(2)/block);
if(numel(rawSize) == 2)
    img1 = LHRforOneChannel(img0, block, blockHeight, blockWidth);
    img1 = uint8(img1);
elseif(numel(rawSize) == 3)
%     img1(:,:,1)= LHRforOneChannel(img0(:,:,1), block, blockHeight, blockWidth); 
%     img1(:,:,2)= LHRforOneChannel(img0(:,:,2), block, blockHeight, blockWidth);
%     img1(:,:,3)= LHRforOneChannel(img0(:,:,3), block, blockHeight, blockWidth);
%     img1 = uint8(img1);
    old_gray_img = rgb2gray(img0);
    new_gray_img = LHRforOneChannel(old_gray_img, block, blockHeight, blockWidth);
    ratio = double(new_gray_img) ./ double(old_gray_img);
    img1(:,:,1)= double(img0(:,:,1)) .* ratio;
    img1(:,:,2)= double(img0(:,:,2)) .* ratio;
    img1(:,:,3)= double(img0(:,:,3)) .* ratio;
    img1 = uint8(img1);
else
    fprintf('wrong img');
end

function[img1] = LHRforOneChannel(img0, block, blockHeight, blockWidth)
cdf = zeros(block * block, 256);
rawSize = size(img0);
for i = 1 : block
    for j = 1 : block
        num = (i-1) * block + j;
        left = (j-1) * blockWidth+1;
        right = j * blockWidth;
        up = (i-1) * blockHeight+1;
        down = i * blockHeight;
        cdf(num,:) = LHRforOneBlock(img0(up:down,left:right));
    end
end
% f1 = cdf(1,:);
% img1(1:uint8(blockHeight/2),1:uint8(blockWidth/2)) = f1(img0(1:uint8(blockHeight/2),1:uint8(blockWidth/2)));
% f2 = cdf(block,:);
% img1(1:uint8(blockHeight/2),uint8((block-1/2)*blockWidth), block*blockWidth) = f2(img0(1:uint8(blockHeight/2),uint8((block-1/2)*blockWidth), block*blockWidth));
% f3 = cdf(block * (block-1) + 1,:);
% img1(uint8((block-1/2)*blockHeight), block*blockHeight,1:uint8(blockWidth/2)) = f3(img0(uint8((block-1/2)*blockHeight), block*blockHeight,1:uint8(blockWidth/2)));
% f4 = cdf(block * block, :);
% img1(uint8((block-1/2)*blockHeight), uint8((block-1/2)*blockWidth), block*blockWidth) = f4(img0(uint8((block-1/2)*blockHeight), block*blockHeight,uint8((block-1/2)*blockWidth), block*blockWidth));

img1 = zeros(rawSize);
for i = 1 : rawSize(1)
    for j = 1 : rawSize(2)
%       four corners
        if i <= blockHeight/2 && j <= blockWidth/2
            f = cdf(1,:);
            img1(i,j) = f(img0(i,j)+1);
        elseif i <= blockHeight/2 && j >= (block-1/2)*blockWidth
            f = cdf(block,:);
            img1(i,j) = f(img0(i,j)+1);
        elseif i >= (block-1/2)*blockHeight && j <= blockWidth/2
            f = cdf(block * (block-1) + 1,:);
            img1(i,j) = f(img0(i,j)+1);
        elseif i >= (block-1/2)*blockHeight && j >= (block-1/2)*blockWidth
            f = cdf(block * block, :);
            img1(i,j) = f(img0(i,j)+1);
%       four edges
        elseif i <= blockHeight/2
            n1 = floor((j - blockWidth/2)/blockWidth) + 1;
            n2 = n1 + 1;
            f1 = cdf(n1,:);
            f2 = cdf(n2,:);
            p =  (j - (n1*blockWidth-blockWidth/2))/blockWidth;  
            q = 1-p;
            img1(i,j) = q * f1(img0(i,j)+1) + p * f2(img0(i,j)+1);
        elseif i >= (block-1/2)*blockHeight
            nY = floor((j - blockWidth/2)/blockWidth);
            n1 = nY + 1 + block * (block-1);
            n2 = n1 + 1;
            f1 = cdf(n1,:);
            f2 = cdf(n2,:);
            p =  (j - ((nY + 1)*blockWidth-blockWidth/2))/blockWidth;  
            q = 1-p;
            img1(i,j) = q * f1(img0(i,j)+1) + p * f2(img0(i,j)+1);
        elseif j <= blockWidth/2
            block_i = floor((i - blockHeight/2)/blockHeight);
            n1 = block_i * block + 1;
            n2 = n1 + block;
            f1 = cdf(n1,:);
            f2 = cdf(n2,:);
            p =  (i - (block_i*blockHeight+blockHeight/2))/blockHeight;  
            q = 1-p;
            img1(i,j) = q * f1(img0(i,j)+1) + p * f2(img0(i,j)+1);
        elseif j >= (block-1/2)*blockWidth
            n1 = block * (floor((i - blockHeight/2)/blockHeight)+1);
            n2 = n1 + block;
            f1 = cdf(n1,:);
            f2 = cdf(n2,:);
            p =  (i - (n1/block*blockHeight-blockHeight/2))/blockHeight;  
            q = 1-p;
            img1(i,j) = q * f1(img0(i,j)+1) + p * f2(img0(i,j)+1);
%       inner space
        else
            block_i = floor((i - blockHeight/2)/blockHeight);
            block_j = floor((j - blockWidth/2)/blockWidth);
            n1 = block_i * block + block_j + 1;
            n2 = n1 + 1;
            n3 = n1 + block;
            n4 = n2 + block;
            p = (j - (block_j*blockWidth+blockWidth/2))/blockWidth;
            q = (i - (block_i*blockHeight+blockHeight/2))/blockHeight; 
            if p < 0 || q < 0 || p > 1 || q > 1
                break;
            end
            f1 = cdf(n1,:);
            f2 = cdf(n2,:);
            f3 = cdf(n3,:);
            f4 = cdf(n4,:);
            img1(i,j) = p*q * f4(img0(i,j)+1) + p*(1-q) * f2(img0(i,j)+1) + (1-p)*q*f3(img0(i,j)+1) + (1-p)*(1-q)*f1(img0(i,j)+1);
        end
    end 
end
img1 = uint8(img1);

function[f] = LHRforOneBlock(inImg)
[pixelCounts, grayLevels] = imhist(inImg);
cdf = cumsum(pixelCounts); % Make transfer function (look up table).
cdf = cdf / sum(pixelCounts); % Normalize
cdfEqual = round(cdf * 255);
alpha = 0.5;
lambda = 1.2;
f = round(alpha * cdfEqual + (1-alpha) * grayLevels);
f = min(f, lambda*grayLevels);