function SFb=compSFb3d(kU,sMaps,opt)
% function to calculate S'*DFT'*U'*kU
% input is the multi-coil kspace
% U1 needs to be multi-coil
% sMaps is same through nt
% 2013/05/07 Yi Guo

% modified 05/20 to 4d


if length(size(kU))>2
    kU_r=kU;
else    
kU_r=zeros(opt.size,opt.class);
kU_r(opt.U)=kU;
end

    b1=conj(repmat(sMaps,[1 1 opt.size(3) 1])).*ifftnd(kU_r,opt.dim);

SFb=sum(b1,4);
end