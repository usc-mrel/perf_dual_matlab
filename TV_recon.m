function imgR = TV_recon(k,Mask,sMaps, lambdas_TV)

%% Parameters for Locally Low Rank (LLR) constraint

lambdas_LLR=[0];  %  LLR set to zero!  turned off!
schatten=0.5;     %[0.6,0.5];
usePar = 0;
overlap = 0;
do_shift = 0;
do_plot = 1;
max_iter = 50;

%% Initalize options from parsed inputs
opt = struct('useParrallel', usePar, 'overlap', overlap, 'do_shift', ...
    do_shift, 'plot', do_plot, 'itermax', max_iter);

%% Sampling Mask
Mask = repmat(Mask,[1 1 1 size(k,1) size(k,5)]);
Mask = permute(Mask,[4 1 2 3 5]);
% sampling pattern U1 the same size as data
U1=(Mask~=0); clear Mask;

[kx,ky,kz,nt,nc]=size(k); 
opt.size=size(k);

opt.B = 4;                      % Patch Size
opt.lambda1=lambdas_TV/10000;	% Spatial(kx,ky) TV penalty
opt.lambda2=lambdas_TV;			% Temporal TV penalty
opt.lambda3=lambdas_LLR;		% LLR penalty
opt.xf=0;  
opt.M=1;                        % spatial mask

if opt.useParrallel
    %Clear any open interactive parrallel computing session
    try
        delete(gcp);
    catch
    end
    parpool(12);            %Setup matlab parrallel pool
end

% build the 3DFT and 3DFT' Fourier operators 
opt.FT = @(im) fft3d (im);
opt.IFT = @(im) ift3d (im);

sMaps=single(sMaps);
sMaps=reshape(sMaps,[kx,ky,kz,1,nc]); %[kx ky kz nc];

%% Recon parameters setting
opt.p1=0.05; 	% dummy variable penaltie from two splitting
opt.p2=0.05;
opt.p3=0.05;
opt.p=schatten;  

opt.class='single';

opt.update=1; % plot update number

opt.U=U1>0;
kU=k(opt.U);

% Use zero-filled images as initial guess
opt.img0=sum(conj(repmat(sMaps,[1 1 1 nt 1])).*opt.IFT(k),5); 

clear k kdata idx header; % clear unnecessary variables to free memory

%% ADMM recon
tic,
[imgR,opt,cost1,cost2]=ADMM2_recon_4d_LLR(kU,sMaps,opt); % run ADMM recon
toc,
imgR=abs(imgR);

end
