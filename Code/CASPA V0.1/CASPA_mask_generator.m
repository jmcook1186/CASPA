

% This function is caled by CASPA_driver.m if the mask functionality is
% switched ON. The idea is for an aerial image of a snowpack to be provided
% in the working directory and the name to be provided when settiung the
% initial condiitons. This script then converts the image to greyscale,
% then binarises and creates a mask, identifying areas of non-snow that are
% 'blocked' in the CASPA cellular model.

function[mask] = CASPA_mask_generator(folder_path, image_name,gridsize)

I = imread(strcat(folder_path,image_name)); %read in each file as an image

B = imresize(I,[sqrt(gridsize), sqrt(gridsize)]);
B = rgb2gray(B);   % Convert indexed to grayscale
level = graythresh(B);  % Compute an appropriate threshold
mask = im2bw(B, level);   % Convert grayscale to binary
mask = not(mask);

save('mask.mat','mask')
