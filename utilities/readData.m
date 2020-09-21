function data = readData(file_name)
% 
%   load data from csv file
%
    opts = detectImportOptions(file_name);
    data = readtable(file_name, opts);
end