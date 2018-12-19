%%%========================Create Radial Bins===========================%%%

% This function gives a cell with the coordinates of element positions for
% the different radial bins. The input arguments are: the size of the image
% in nm; the size of the matrix in pixels (this is input as any matrix of
% the correct size); the radius for the first bin, in nm; the radius for
% all subsequent bins; and finally, the number of bins in total. The
% outputs are: the coordinates of the bins as a cell; the matrix with the
% nm value corresponding to each pixel; and an array containing the average
% radial value, in nm, for each bin.

function [Circle_coords, Circle_Matrix, Radial_Value] = RadialBins(Crop_Size, Matrix_Of_Same_Size, Initial_Radius, Radius, NumberOfBins)

% Rows and columns of a cropped image
[rr, cc] = size(Matrix_Of_Same_Size);

Ones_Matrix = ones(rr,cc);

% nm/pixels for 130nm crop
Ratio_nm_per_pixs = Crop_Size/rr;

% Pre-allocate cell array for the coordinates of the pixels of the circles.
% The number of cells will be the number of bins.
Circle_coords      = cell(1,NumberOfBins);
Mask_Circle        = cell(1,NumberOfBins);
Mask_Store         = cell(1,NumberOfBins);  
Mask_Circle_Pixels = cell(1,NumberOfBins);
Radial_Value       = zeros(1,NumberOfBins);

% Manually create the first bin (0-6nm radius)
min_boundary = 0;
max_boundary = Initial_Radius;
[Mesh_x, Mesh_y] = meshgrid(-((rr/2)-0.5):((cc/2)-0.5));
% C is a matrix whose central element = 0, and all elements leading away
% increase linearly in a circular pattern.
C = sqrt(((Mesh_x).^2) + ((Mesh_y).^2));
% Circle matrix converts these values using the scalar Ratio_nm_per_pixels
% to convert the matrix to reflect the real nm values.
Circle_Matrix = C*Ratio_nm_per_pixs;
% Index the nm values desired.
Mask_Circle{1} = Circle_Matrix>=min_boundary & Circle_Matrix<max_boundary;
Mask_Circle_Pixels{1} = Circle_Matrix(Circle_Matrix >=...
        min_boundary & Circle_Matrix < max_boundary);
% Use the mask to delete all other pixel values.
Masked_Pore = Ones_Matrix.*Mask_Circle{1};

[Coords_rr, Coords_cc] = find(Masked_Pore>0);
Circle_coords{1} = [Coords_rr, Coords_cc];


Mask_Store{1} = zeros(rr,cc);
for j=1:length(Coords_rr)
    Mask_Store{1}(Coords_rr(j),Coords_cc(j))=1;
end

for i=1:NumberOfBins-1
    Min_Boundary = Initial_Radius + ((i-1)*Radius);
    Max_Boundary = Initial_Radius + (i*Radius);
    Mask_Circle{i+1} = Circle_Matrix >=...
        Min_Boundary & Circle_Matrix < Max_Boundary;
    Mask_Circle_Pixels{i+1} = Circle_Matrix(Circle_Matrix >=...
        Min_Boundary & Circle_Matrix < Max_Boundary);
    Masked_Ones = Ones_Matrix.*Mask_Circle{i+1};
    [Coords_rr, Coords_cc] = find(Masked_Ones>0);
    Circle_coords{i+1} = [Coords_rr, Coords_cc];
    Mask_Store{i+1} = zeros(rr,cc);
    for j=1:length(Circle_coords{i+1})
        Mask_Store{i+1}(Coords_rr(j),Coords_cc(j))=1;
    end
end


% Need to save values of pixels in each bin, radial values of circle
% matrix, so can find average radial value for the radii array.
for i=1:NumberOfBins
     Radial_Value(i) = mean(Mask_Circle_Pixels{i});
end

% overwrite first value and set to 0nm.
Radial_Value(1) = 0;

%%%=====================End Create Radial Bins===========================%%%

end