function res=compSx3d(x,sMaps,opt)
res=repmat(sMaps,[1 1 opt.size(3) 1]).*repmat(x,[1 1 1 opt.size(4)]);
end