function [kd_new,idxd_new, ks_new,idxs_new] = viewsharing(kd,idxd,ks,idxs,step)

[nx ny nz nt nc] = size(kd);
nt_new = (nt-step+1);
kd_new = zeros(nx,ny,nz,nt_new,nc);
idxd_new = zeros(ny,nz,nt_new);

[nx ny nz nt nc] = size(ks);
ks_new = zeros(nx,ny,nz,nt_new,nc);
idxs_new = zeros(ny,nz,nt_new);

for nn=1:nt_new
    
    kd_new(:,:,:,nn,:) = sum(kd(:,:,:,nn:(nn+step-1),:),4);
    idxd_new(:,:,nn) = sum(idxd(:,:,nn:(nn+step-1)),3);
    
    ks_new(:,:,:,nn,:) = sum(ks(:,:,:,nn:(nn+step-1),:),4);
    idxs_new(:,:,nn) = sum(idxs(:,:,nn:(nn+step-1)),3);
end

% tmp_d = repmat(idxd_new,[1 1 1 nx nc]);
% tmp_d = permute(tmp_d,[4 1 2 3 5]);
% kd_new(tmp_d~=0) = kd_new(tmp_d~=0)./tmp_d(tmp_d~=0) ;
% 
% tmp_s = repmat(idxs_new,[1 1 1 nx nc]);
% tmp_s = permute(tmp_s,[4 1 2 3 5]);
% ks_new(tmp_s~=0) = ks_new(tmp_s~=0)./tmp_s(tmp_s~=0) ;
% 
% idxd_new = logical(idxd_new);
% idxs_new = logical(idxs_new);

end


