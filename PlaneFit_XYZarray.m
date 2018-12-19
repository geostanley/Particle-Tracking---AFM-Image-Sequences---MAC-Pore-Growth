%%%=== PlaneFit_XYZarray ===%%%

% This function takes in an Nx3 array (i.e., XYZ coordinates), and fits a
% plane to the data (output as a matrix the same size as 'matrix', from
% which the XYZarray would have been created).


function [plane] = PlaneFit_XYZarray(matrix, XYZ_array)

    X_array = XYZ_array(:,1);
    Y_array = XYZ_array(:,2);
    Z_array = XYZ_array(:,3);

    C = planefit(X_array, Y_array, Z_array);

    a = C(1);
    b = C(2);
    c = C(3);

    plane = zeros(size(matrix));

    [row,col]=size(matrix);
    
    for i=1:col
        for j=1:row
            plane(j,i) = (i*a) + (j*b) + c;
        end
    end
    
end