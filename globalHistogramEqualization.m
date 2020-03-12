function[img1] = globalHistogramEqualization(img0)
%Global histogram equalization for 2D colorful image
rawSize = size(img0);
if(numel(rawSize) == 2)
    img1 = GHRforOneChannel(img0);
elseif(numel(rawSize) == 3)
%     img1(:,:,1)= GHRforOneChannel(img0(:,:,1)); 
%     img1(:,:,2)= GHRforOneChannel(img0(:,:,2));
%     img1(:,:,3)= GHRforOneChannel(img0(:,:,3));
%     img1 = uint8(img1);
    oldGrayImg = rgb2gray(img0);
    newGrayImg = GHRforOneChannel(oldGrayImg);
    ratio = double(newGrayImg) ./ double(oldGrayImg);
    img1(:,:,1)= double(img0(:,:,1)) .* ratio;
    img1(:,:,2)= double(img0(:,:,2)) .* ratio;
    img1(:,:,3)= double(img0(:,:,3)) .* ratio;
    img1 = uint8(img1);
else
    fprintf('wrong img');
end

function[img1] = GHRforOneChannel(img0)
%globalHistogramEqualization for 2D image
%img1: Equalization picture
%img0: raw picture
[pixelCounts, grayLevels] = imhist(img0);
cdf = cumsum(pixelCounts); % Make transfer function (look up table).
cdf = cdf / sum(pixelCounts); % Normalize
cdfEqual = round(cdf * 255);
alpha = 0.9;
f = round(alpha * cdfEqual + (1-alpha) * grayLevels);
lambda = 1.2;
f = min(f, lambda*grayLevels);
picEqual = f( img0 + 1 );
img1 = uint8(picEqual);