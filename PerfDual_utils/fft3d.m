function f = fft3d(f,shift)
%   Computes the [shifted] 3D fft

if nargin < 2 || isempty(shift)
    shift = 1;
end

if shift
    f = ifftshift(ifftshift(ifftshift(fft(fft(fft(fftshift(fftshift(fftshift(...
        f,1),2),3),[],1),[],2),[],3),1),2),3);
%     f = fftshift(fftshift(fftshift(fft(fft(fft(f,[],1),[],2),[],3),1),2),3);
else
    f = fft(fft(fft(f,[],1),[],2),[],3);
end
