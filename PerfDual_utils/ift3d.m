function F = ift3d(F,shift)
%   Computes the [shifted] 3D ifft

if nargin < 2 || isempty(shift)
    shift = 1;
end

if shift
    F = fftshift(fftshift(fftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(...
        F,1),2),3),[],1),[],2),[],3),1),2),3);
%     F = ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(F,1),2),3),[],1),[],2),[],3);
else
    F = ifft(ifft(ifft(F,[],1),[],2),[],3);
end
