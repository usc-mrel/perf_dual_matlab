function Shx=compShx3d(Is,sMaps,opt)
% function to calculate Sh*x
% input Is is the multi-coil images
% U1 needs to be multi-coil
% sMaps is same through nt

% 2013/05/07 Yi Guo


x=conj(repmat(sMaps,[1 1 opt.size(3) 1])).*Is;
Shx=sum(x,4);
end