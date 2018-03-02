function [ mask sampidx] = mask_lattice( xx,yy,zz,accel_ky, accel_kz )
%function [ mask sampidx] = mask_lattice( xx,yy,zz,accel_ky, accel_kz )
%
% Input
%    xx - freq encoding dimension
%    yy - phase encoding dimension
%    zz - partition encoding dimension
%    accel_ky - acceleration factor in phase encoding dimension
%    accel_kz - acceleration factor in partition encoding dimension
%
% Output
%    mask - binary mask matrix [xx,yy,zz, accel]
%    sampidx - indice of sampled position based on 1D vector format of 3D
%    data (xx,yy,zz)
%
%              Taehoon Shin (shinage@gmail.com)  Mar,2008

if (nargin < 6 )
    undtype = 0;
end


if ( mod(yy,accel_ky) ~= 0 )
    error('accel_ky should be divisor of yy');
end

if ( mod(zz,accel_kz) ~= 0 )
    error('accel_kz should be divisor of zz');
end

accel = accel_ky * accel_kz;

% Constructing 3D binary mask
mask = zeros(xx,yy,zz, accel);
idx = 0;
for ii=1: accel_ky
    for jj=1:accel_kz
        idx = idx + 1;
        mask(:,ii:accel_ky:end, jj:accel_kz:end, idx )=1;
    end
end

% Store indice for sampled position in 2D ky-kz plane
for kk=1: accel
    idx = 0;
    for aa=1: yy
        for bb=1:zz
            if mask(1, aa,bb, kk) == true
                idx = idx + 1;
                sampidx(idx, kk) =  aa + i*bb;
            end
        end
    end
end

            
    
    
   