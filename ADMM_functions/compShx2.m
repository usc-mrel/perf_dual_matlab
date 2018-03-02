function img = compShx2(imgS,sMaps,size)
%   Applies the inverse SENSE operator to the image
%   
%   Author: RML
%   Date: 11/2011
%   
%   Usage: img = iSENSE(imgS,opt)
%   
%   Input:
%   imgS: Coil images of size RO x PE x NS x NT x NR
%   opt: Specifies transform options for recon. Use SPIRiT_optset.m
%   	Must include fields:
%           opt.S: Sensitivity maps of size:
%               RO x PE x NS x 1 x NR
%   
%   Output:
%   img: image of size RO x PE x NS x NT

% %   Get number of time frames
% nt = opt.size(4);
% 
% %   Combine coils into singe image
% if isfield(opt,'Si') && ~isempty(opt.Si)
%     img = sum(imgS .* repmat(opt.Si,[1 1 1 nt 1]),5);
% else
%     img = sum(imgS .* repmat(conj(opt.S),[1 1 1 nt 1]),5);
% end

img = iSENSE_MEX(imgS,sMaps,size);

end
