function [Rx,fsize,LL]=compRx4d(x,opt)
% Overall sparse transform 4d wavelet + 3d TV + Temporal filter
% add automatic wavelet size organizer
% add dimension and constraints detector

% Yi Guo 
% 07/31/2013

% x=reshape(x,opt.kx,opt.ky,opt.kz,opt.nt);

if opt.beta1~=0

wc=fwtN(x,'db2',opt.worder); % get wavelet coefficient

Wx=[]; % stack all coefficient to this Wx vector
N = length(wc);
fsize=cell(1,N);

%   Loop through transform order and the various coefficients
for i = 1:N
    WCt = wc{i};
    fnames = fieldnames(WCt);
    nf = size(fnames);
    for j = 1:nf
        fname = fnames{j};
        if ~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ~strcmp(fname,'nd')
            coefs = WCt.(fname);
            if ~isempty(coefs)
                Wx=cat(1,Wx,coefs(:));
            end
            fsize{i}.(fname)=size(coefs); % store wavelet name and size in fsize
            
        else
            fsize{i}.(fname)=WCt.(fname);
        end
    end
end

LL(1)=length(Wx);

else
LL(1)=0;
Wx=[];
fsize=0;
end

if opt.beta2~=0
    if opt.size(3)>1

% Calculate Total variation using shift difference
TV1=x-circshift(x,[1 0 0 0]);
TV2=x-circshift(x,[0 1 0 0]);
TV3=x-circshift(x,[0 0 1 0]);
TV=[TV1(:);TV2(:);TV3(:)];
LL(2)=numel(x)*3;
    end
    
    if opt.size(3)==1
TV1=x-circshift(x,[1 0 0 0]);
TV2=x-circshift(x,[0 1 0 0]);
TV=[TV1(:);TV2(:)];
LL(2)=numel(x)*2;
    end
else
    LL(2)=0;
    TV=[];
end

if opt.beta3~=0 && opt.size(4)>5
% linear high-pass filter with Marc's kernel: [0.06 0.44 -1 0.44 0.06]
% Hx=0.06*x+0.44*circshift(x,[0 0 0 1])-circshift(x,[0 0 0 2])+0.44*circshift(x,[0 0 0 3])+0.06*circshift(x,[0 0 0 4]);

% TV in the temporal dimension
Hx=x-circshift(x,[0 0 0 1]);

% Hx(:,:,:,1)=0;
%[kx,ky,kz,nt]=size(x);
%W=cat(4,zeros(kx,ky,kz,1),ones(kx,ky,kz,5)*1,ones(kx,ky,kz,3)*0.1,ones(kx,ky,kz,26)*2);
%Hx=W.*Hx;
LL(3)=numel(Hx);
else
    LL(3)=0;
    Hx=[];
end

Rx=[Wx;TV;Hx(:)];

end