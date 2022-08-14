function testDiffusion ()
num_iter = 6;
delta_t = 0.75;
kappa = 2;
option = 1;


im1 = rgb2gray(imread('brainweb59.tif'));
%im1 = (dicomread('slice_3D.dcm'));
Im = imgaussfilt(im1,5);
Im = im2double(Im)
%montage({Im,im1})
%title('Original Image (Left) Vs. Gaussian Filtered Image (Right)')

%blurring
H = fspecial('disk',3);
blurred = imfilter(im1,H,'replicate');

ad = anisodiff(Im,num_iter,kappa,delta_t,option);
figure, subplot 131, imshow(im1), subplot 132, imshow(ad,[]), subplot 133, imshow(blurred,[])



