%% INITIALIZATION

clear all; close all;

% adding local subdirectories to the path (supporting functions)
addpath('SKED');
addpath('ADMM_functions');
addpath('PerfDual_utils');
addpath('sample_data');

%% IMAGING PARAMETERS
ncoil = 8;    % number of receiver coils
accel_ky = 3; % undersampling factor along k_y
accel_kz = 2; % undersampling factor along k_z

flag_viewshare = 1;  % perform view-sharing (borrow data from temporally adjacent time points)
lambda = 0.1;        % regularization parameter (how much to weight "total variation sparsity")

%% LOAD DATA AND SETUP DATA STRUCTURE

[dr,rhfrsize,usercv,rhnframes,rhnslices,hdr] = rawloadHD_si('P39936.7');  % loads GE Raw Data File (a.k.a. "P-file")

[kd,idxd,ks,idxs,usercv] = sampling(dr,rhfrsize,usercv,rhnframes,rhnslices,ncoil,accel_ky,accel_kz); % generates the ky,kz,t sample locations

if flag_viewshare == 1 % For quick look (no reconstruction constraints, view sharing)
    [kd_vs, ~, ks_vs, ~] = viewsharing(kd, idxd, ks, idxs, accel_ky*accel_kz);
    
    % Diastole
    imgR_D = ift3d(kd_vs);
    imgR_D = sum(imgR_D.*conj(imgR_D), 5);
    imgR_D = imresize_perf(imgR_D, usercv, 1); 
    im1 = cat(2, imgR_D(:,:,1,:), imgR_D(:,:,2,:), imgR_D(:,:,3,:), imgR_D(:,:,4,:), imgR_D(:,:,5,:));
    im2 = cat(2, imgR_D(:,:,6,:), imgR_D(:,:,7,:), imgR_D(:,:,8,:), imgR_D(:,:,9,:), imgR_D(:,:,10,:));
    im = cat(1,im1,im2);
    im = im/3;
    im(im>1)=1;
    img2mp4( abs(im), 10,'Diastole_viewsharing.mp4');
    
    % Systole
    imgR_S = ift3d(ks_vs);
    imgR_S = sum(imgR_S.*conj(imgR_S), 5);
    imgR_S = imresize_perf(imgR_S,usercv,0);
    im1s = cat(2, imgR_S(:,:,1,:), imgR_S(:,:,2,:), imgR_S(:,:,3,:), imgR_S(:,:,4,:));
    im2s = cat(2, imgR_S(:,:,5,:), imgR_S(:,:,6,:), imgR_S(:,:,7,:), imgR_S(:,:,8,:));
    ims = cat(1,im1s,im2s);
    ims = ims/10;
    ims(ims>1)=1;
    img2mp4( abs(ims), 10,'Systole_viewsharing.mp4');
end

clear imgR_D imgR_S im ims;

%% ADMM for diastolic acq data
[sMaps] = sMaps_sos_4Ddata(kd, idxd);

imgR_D = TV_recon(kd,idxd,sMaps,lambda);

imgR_D = imresize_perf(imgR_D,usercv,1); %diastole

im1 = cat(2, imgR_D(:,:,1,:), imgR_D(:,:,2,:), imgR_D(:,:,3,:), imgR_D(:,:,4,:), imgR_D(:,:,5,:));
im2 = cat(2, imgR_D(:,:,6,:), imgR_D(:,:,7,:), imgR_D(:,:,8,:), imgR_D(:,:,9,:), imgR_D(:,:,10,:));
im = cat(1,im1,im2);

[nx,ny,nz,nt]=size(imgR_D);
figure;montage(reshape(abs(imgR_D(:,:,:,15)),[nx,ny,1,nz]),'DisplayRange',[0 3])

im = im/3;
im(im>1)=1;
img2mp4( abs(im), 10,'Diastole_temporal_TV.mp4');
%% ADMM for systolic acq data

[sMaps] = sMaps_sos_4Ddata(ks, idxs);

imgR_S = TV_recon(ks,idxs,sMaps,lambda);

imgR_S = imresize_perf(imgR_S,usercv,0);%systole

im1s = cat(2, imgR_S(:,:,1,:), imgR_S(:,:,2,:), imgR_S(:,:,3,:), imgR_S(:,:,4,:));
im2s = cat(2, imgR_S(:,:,5,:), imgR_S(:,:,6,:), imgR_S(:,:,7,:), imgR_S(:,:,8,:));
ims = cat(1,im1s,im2s);

[nx,ny,nz,nt]=size(imgR_S);
figure;montage(reshape(abs(imgR_S(:,:,:,15)),[nx,ny,1,nz]),'DisplayRange',[0 10])

ims = ims/10;
ims(ims>1)=1;
img2mp4( abs(ims), 10,'Systole_temporal_TV.mp4');


    

