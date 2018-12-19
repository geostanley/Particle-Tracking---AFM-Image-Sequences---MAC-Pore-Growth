% This function uses the ginput function to allow the user to select a
% template. This can then be used with the normxcorr2 function to pick out
% similar features in an image.

% input the entire image file (the one from which you want to select a
% template); the width of the cropped template, in nm; and, the nm/pixel
% value

% the output is a matrix of your chosen template. The width of the template
% is adjusted to ensure an odd number of pixels. This ensures a
% well-defined centre of the template (for rotational averaging - if
% required).

function [template] = SelectTemplate_even(heightdata, CropSize_nm, nmperpixel, clims)

% select template
centres = InputCentres_Copper_CX_manual_clims(heightdata, clims);

round_centres = round(centres);

% get the coordinates for the corner of the template to be cropped
Crop_Size_Half = CropSize_nm/2;
c1 = round_centres(1,1)-Crop_Size_Half/nmperpixel;
c2 = round_centres(1,2)-Crop_Size_Half/nmperpixel;

% ensure an even number of pixels along edges
numb_pixels_crop = round(CropSize_nm/nmperpixel);
if mod(numb_pixels_crop,2) == 0 % if numb of pixels is odd, remove one to make even
    numb_pixels_crop = numb_pixels_crop - 1;
end

% size1 = size(imcrop(heightdata, [c1 c2 numb_pixels_crop numb_pixels_crop]),1);
% size2 = size(imcrop(heightdata, [c1 c2 numb_pixels_crop numb_pixels_crop]),2);


template = imcrop(heightdata, [c1 c2 numb_pixels_crop numb_pixels_crop]);


end

