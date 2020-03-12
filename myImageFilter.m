function[img1] = myImageFilter(img0, k, s)
%convolution for 2D image
%img1: convolution picture
%img0: raw picture
%k: kernel
%s: stride, e.g.[1,1], which means shifting step by step both vertically and horizonally
if nargin < 2
    k = [1/9 1/9 1/9;1/9 1/9 1/9;1/9 1/9 1/9];
    s = [1 1];
elseif nargin < 3
    s = [1 1];
end
[h, w, c] = size(img0);
[kerh, kerw] = size(k);
paddingImg = zeros(h + 2*(kerh - 1), w + 2*(kerw - 1), c);
paddingImg(kerh - 1: h + kerh - 2, kerw - 1: w + kerw - 2, :) = img0;
convSize = [ceil((h + kerh -1)/s(1)),ceil((w + kerw - 1)/s(2))];
convImg = zeros([convSize, c]);
k = repmat(k, [1,1,c]);
rowIdx = 1 : s(1) : size(paddingImg, 1);
colIdx = 1 : s(2) : size(paddingImg, 2);
for row = 1 : convSize(1)
    up = rowIdx(row);
    down = up + kerh - 1;
    for col = 1 : convSize(2)
        left = colIdx(col);
        right = left + kerw - 1;
        convImg(row,col,:) = sum(sum(paddingImg(up:down,left:right,:).*k,1),2);
    end
end
img1 = uint8(convImg);
end