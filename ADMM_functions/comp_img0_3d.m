function img=comp_img0_3d(k,sMaps,opt)
% calculate time-averaged image of 3d data (kx,ky,nt)

kU = zeros(opt.size,'single');
kU(opt.U) = k;
km = mean(kU,3)./(mean(kU~=0,3) + eps);

km= repmat(km,[1 1 opt.size(3) 1]);
imgM = opt.IFT(km);

img=compShx3d(imgM,sMaps,opt);

end