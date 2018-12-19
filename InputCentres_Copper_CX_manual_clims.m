% This function is used to assign the coordinates for the central axes of
% rotation for the NPCs. The input argument is the AFM image in greyscale.
% A figure is then produced showing the image, and the user must manually
% click on the picture to assign the central axes of rotation for each NPC.
% The output is an Nx2 array of the coordinates, to the nearest pixel, of
% the selected centres.

function Centres = InputCentres_Copper_CX_manual_clims(Image, clims)

hFig = figure(); imagesc(Image, clims);
set(hFig, 'Position', [300 10 1000 1000])

colormap(copper);
xlabel('Select centre and press enter when finished')
title('Select a template MAC and press enter when finished', 'FontSize', 14)
set (gcf, 'WindowButtonMotionFcn', @mouseMove);
[c, r] = ginput;
r = round(r);
c = round(c);
Centres = [c r];
       
    
end