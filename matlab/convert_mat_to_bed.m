data_dir = '../data/k/';
mode = 'w';
chrs = [1:3 10:22];
for i = 1: 16
    data_fp = strcat(data_dir, 'MLRates_Protocol1_Chr', int2str(chrs(i)),'WT.mat');
    data = load(data_fp);
    k = data.MLELam(:, 1);
    sites = data.FitSites(:,1);
    
    if i > 1
        mode = 'a';
    end
    matrix = ones(length(sites), 4)*chrs(i);
    matrix(:, 2) = sites;
    matrix(:, 3) = sites;
    matrix(:, 4) = k;
    
    fileID = fopen('../data/K.bed', mode);
    fprintf(fileID,'chr%d\t%d\t%d\t%.4f\n', matrix');
    fclose(fileID);
end