function ind = patchInd(X,B,overlap,varargin)
%test for patch blocks
% Remember to add in varargin (for mask to limit overlapping patches)

%Obtain Dimensions
%[nx, ny, n3, n4] = size(X);
[nx, ny, nz, n3] = size(X);
nx_mod = nx-mod(nx,B);
ny_mod = ny-mod(ny,B);
nz_mod = nz-mod(nz,B);
%Overlapping patches is performed by simply applying the amount of 
%shift to the x and y indices

if overlap == 1
    [shiftx, shifty] = meshgrid(linspace(0,B-1,3), linspace(0,B-1,3));
    shiftx = round(shiftx(:));
    shifty = round(shifty(:));
else
    shiftx = 0; shifty = 0; shiftz = 0;
end

Nshift = length(shiftx);

%Find linear indices within the nx and ny dimension 

%Patch indices in x-dim
x = repmat([1:B]',[B 1]);
x = repmat(x, [ny_mod/B 1]);
inc = B*ones(numel(x), floor(nx_mod/B)-1);
x = [x inc];
x = cumsum(x,2);
x=x(:);

%Patch indicies in y-dim
y = ones(B, ny_mod);
y = cumsum(y,2);
y = repmat(y(:), [floor(nx_mod/B) 1]);

z = ones(B, nz_mod);
z = cumsum(z,2);
z = repmat(z(:), [floor(nx_mod/B) 1]);

ind = zeros(nx_mod*ny_mod*nz_mod*n3,Nshift);

% Stack all the x and y indices including shift for
% overlapping patches
for k = 1:Nshift
	% Convert to linear indicies 
	sub_ind = sub2ind([nx ny nz], mod(x+shiftx(k),nx)+1, mod(y+shifty(k),ny)+1,mod(z+shiftz(k),nz)+1);

	%Replicate into the third and fourth dimension
	%ind = bsxfun(@plus,sub_ind,(0:(n3-1)*(n4-1)).*(nx*ny));

	%Replicate 2D linear indices to 3rd dimesion
	temp_ind = bsxfun(@plus,sub_ind,(0:(n3-1)).*(nx*ny*nz));

	%Repliccate 3D linear indcies into 4th dimension
	%ind = bsxfun(@plus,ind(:), (0:(n4-1)).*(nx*ny*n3));
	ind(:,k) = temp_ind(:);
end
