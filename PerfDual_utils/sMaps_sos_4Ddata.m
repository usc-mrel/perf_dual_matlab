function [sMaps] = sMaps_sos_4Ddata (k, Mask)

k = permute(k,[1 2 3 5 4]);

k = sum(k,5); %average along temporal dimension

Mask = repmat(sum(Mask,3),[ 1 1 size(k,1) size(k,4)]);
Mask = permute(Mask,[3 1 2 4]);

k(Mask~=0) = k(Mask~=0)./Mask((Mask~=0));
Kdata = k; clear k;
%% ifft in the readout direction, prepare to smaps slice by slice
[nreadout,ky,kz,ncoil]=size(Kdata);

Im=ift3d(Kdata); clear Kdata;

Im_sum = sqrt(sum(Im.*conj(Im),4));
for nnc = 1:ncoil
    sMaps(:,:,:,nnc) = Im(:,:,:,nnc)./Im_sum;
end

tmp=reshape(abs(sMaps),[nreadout,ky,kz*ncoil]);
figure;montage(reshape(tmp,[nreadout,ky,1, kz*ncoil]),'DisplayRange',[]);

end
