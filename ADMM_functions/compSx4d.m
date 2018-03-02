function Sx=compSx4d(x,sMaps,opt)
% function to calculate S*x
% input is the single frame(coil) image
% output is the multi frame(coil) images
% 2013/05/07 Yi Guo

% update 05/17 to do 3d
% 05/20 to 4d

% x_k=reshape(x,opt.kx,opt.ky,opt.kz,opt.nt);
% Sx=zeros(opt.kx,opt.ky,opt.kz,opt.nt,opt.coils);


Sx=repmat(sMaps,[1 1 1 opt.size(4) 1]).*repmat(x,[1 1 1 1 opt.size(5)]);

end
