function [tau, nb, na, K, max_iter, error, alghorithm, learning_type] = readConfig()
%readConfig 
%   Detailed explanation goes here
    %Save to file
    file_template = '%i %i %i %i %i %f %i %i';
    fileID = fopen('ustawienia.txt','r');
    temp = fscanf(fileID,file_template);
    fclose(fileID);
    
    tau = temp(1);
    nb = temp(2);
    na = temp(3);
    K = temp(4);
    max_iter = temp(5);
    error = temp(6);
    alghorithm = temp(7);
    learning_type = temp(8);
end

