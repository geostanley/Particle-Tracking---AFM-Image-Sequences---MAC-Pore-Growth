%%%=== Create_FileNames_Cell_2 ===%%%

% This function creates a cell containing strings of filenames. Input the
% directory contatining the desired data files, the file suffix common to
% all data files, and an array containing the numbers of all file ends.

% Also tell the function which version of Nanoscope the spm files are from.
% If 9.1, input the string 'yes'. If 9.2 or higher, input 'no'.

% A cell array with all the file names will be returned.

function [FileNames] = Create_FileNames_Cell_2(Data_Directory, Generic_FileName, File_Nos, NS_nine_one)

File_Numbers = cell(1,length(File_Nos));
for i=1:length(File_Nos);
    if File_Nos(i) < 10
        File_Numbers{i} = strcat('00', num2str(File_Nos(i)));
    elseif File_Nos(i) >=10 && File_Nos(i) < 100
        File_Numbers{i} = strcat('0', num2str(File_Nos(i)));
    else
        File_Numbers{i} = num2str(File_Nos(i));
    end
end


if length(NS_nine_one) == 3

    FileNames = cell(1,length(File_Nos));

    for i = 1:length(File_Nos)

        fullfilename = fullfile(Data_Directory, Generic_FileName);
        FileNames{i} = horzcat(fullfilename,'.', File_Numbers{i});

    end
    
else
    
    FileNames = cell(1,length(File_Nos));

    for i = 1:length(File_Nos)

        fullfilename = fullfile(Data_Directory, Generic_FileName);
        FileNames{i} = horzcat(fullfilename,'.0_00', File_Numbers{i}, 'Flatten Image View_1.spm');

    end
    
end

end