function Shx=compShx4d(Is,sMaps,opt)
% function to calculate Sh*x
% input Is is the multi-coil images
% U1 needs to be multi-coil
% sMaps is same through nt

% 2013/05/07 Yi Guo

% modified to do 4d 05/20

% Is=reshape(Is,opt.kx,opt.ky,opt.kz,opt.nt,opt.coils);

x=conj(repmat(sMaps,[1 1 1 opt.size(4) 1])).*Is;

Shx=sum(x,5);
% Shx=x(:);
end