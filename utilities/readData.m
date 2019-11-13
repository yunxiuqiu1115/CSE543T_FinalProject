function data = readData(file_name)
    opts = detectImportOptions(file_name);
    data = readtable(file_name, opts);
end