function [skedky,skedkz] = readSkedFile(skedfile)


    fid = fopen(skedfile, 'r');

    %% read sked file
    opfov = fscanf(fid, '%d', 1);
    opzfov = fscanf(fid, '%d', 1);
    opxres = fscanf(fid, '%d', 1);
    opyres = fscanf(fid, '%d', 1);
    opzres = fscanf(fid, '%d', 1);
    nsegs = fscanf(fid, '%d', 1); % total number of views
    nsegs_accel = fscanf(fid, '%d', 1);% nsegs_accel = nsegs/R

    skedky = zeros(nsegs,1);
    skedkz = zeros(nsegs,1);

    for i = 1: nsegs
        skedky(i)=fscanf(fid, '%d', 1);
        skedkz(i)=fscanf(fid, '%d', 1);
    end

end