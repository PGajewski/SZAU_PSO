

function [u,y] = readData(file)
%READDATA Summary of this function goes here
%   Detailed explanation goes here
    file_template = '%f %f\n';
    fileID = fopen(file,'r');
    temp = fscanf(fileID,file_template, [2,inf]);
    fclose(fileID);
    
    u = temp(1,:);
    y = temp(2,:);
end

