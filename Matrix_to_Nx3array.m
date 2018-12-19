%%%=== Matrix_to_Nx3array ===%%%

% This function takes in a matrix (square or rectangular), and converts 
% it into an Nx3 array. I.e., XYZ coordinates.

function [XYZ_array] = Matrix_to_Nx3array(matrix)

    % get the dimensions of the cropped image. As square, r=c.
    [r,c] = size(matrix);

    % transform the square matrix into an Nx3 array
    XYZ_array = zeros(length(matrix(:)), 3);

    for i=1:c
        
        XYZ_array((((i-1)*r+1):i*r), 1) = i;
        
    end
    
    
    for j = 1:c
        
        XYZ_array((((j-1)*r+1):j*r), 2) = 1:r;
            
    end

    for i=1:length(XYZ_array)
        XYZ_array(i,3) = matrix(XYZ_array(i,2),XYZ_array(i,1));
    end
    
    
end