function imgS = compSx2(img,sMaps,size)
%   Applies the forward SENSE operator to the image
%   
%   Author: RML
%   Date: 11/2011
%   
%   Usage: imgS = fSENSE(img,opt)
%   
%   Input:
%   img: Combined coil image of size RO x PE x NS x NT x 1
%   opt: Specifies transform options for recon. Use SPIRiT_optset.m
%   	Must include fields:
%           opt.S: Sensitivity maps of size:
%               RO x PE x NS x 1 x NR
%   
%   Output:
%   img: image of size RO x PE x NS x NT x NR

% %   Get sizes
% nt = opt.size(4);
% nr = opt.size(5);

%   Make multicoil image
% imgS = repmat(img,[1 1 1 1 nr]).*repmat(opt.S,[1 1 1 nt 1]);

% %   Alternate approach (no repmat)
% imgS = zeros(opt.size,opt.class);
% for i = 1:nt
%     for j = 1:nr
%         imgS(:,:,:,i,j) = img(:,:,:,i).*opt.S(:,:,:,1,j);
%     end
% end

imgS = fSENSE_MEX(img,sMaps,size);

end
