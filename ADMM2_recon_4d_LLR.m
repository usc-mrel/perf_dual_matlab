function [imgR,opt,cost1,cost2,cost3]=ADMM2_recon_4d_LLR(k,sMaps, opt)
% function to realize 4d CS SENSE using ADMM which has two splitting steps
% input: kU: undersampled kspace (zero-filled)
%        sMaps: sensitivity Map, nt dimension should be 1, even for 4d
%        U1: sampling pattern, same size with kU
%        opt: recon parameters
% output: imgR: reconstructed image(kx*ky*kz)
%         opt: other necessary output

% Yi Guo 07/29/2013

kx=opt.size(1);ky=opt.size(2);kz=opt.size(3);nt=opt.size(4);coils=opt.size(5);
% [kx,ky,kz,nt,coils]=size(kU); % get data size

%% set initial values
x_k=opt.img0;
%im=squeeze(abs(x_k(4,:,:,:)));im=double(im);ims2avi('initial.avi',0.1,im/max(im(:)),0);
% x_kp=ones(kx*ky*kz*nt,1,'single');

[v1_k,LL]=compTx4d(x_k,opt); % v=Rz

v2_k=x_k;
u_k=compSx4d(x_k,sMaps,opt);

e11_k=zeros(LL,1,opt.class); %residue variable for TV
e12_k=zeros(kx,ky,kz,nt,opt.class); % residue variable for LLR
e2_k=zeros(kx,ky,kz,nt,coils,opt.class); % size of k-space
e3_k=zeros(kx,ky,kz,nt,opt.class); % size of image

% set different weighting for spatial and temporal TV
if ((opt.lambda1==0) &&  (opt.lambda2~=0))
    TVlambda=[ones([kx*ky*kz*nt,1],opt.class)*opt.lambda2];
    disp('temporal TV');
elseif ((opt.lambda1~=0) &&  (opt.lambda2==0))
    TVlambda=[ones([kx*ky*kz*nt*3,1],opt.class)*opt.lambda1];
    disp('spatial TV');
else
    TVlambda=[ones([kx*ky*kz*nt*3,1],opt.class)*opt.lambda1;ones([kx*ky*kz*nt,1],opt.class)*opt.lambda2];
    disp('temporal+spatial TV');
end

%% Pre-calculate the inverse matrix
SIinv=opt.p2*repmat(sum(conj(sMaps).*sMaps,5),[1 1 1 nt])+opt.p3;
UIinv=opt.U+opt.p2;
Tinv=invQT4d(opt);

%%
iter=1;

tic,
tc=cputime;

cost=zeros(1,opt.itermax);

while iter<opt.itermax
tic

x_kp=x_k; % store previous x_k

%---------------------------------------------------------------------------
% if opt.lambda1~=0 && opt.lambda3~=0
% z_k=ifftn(fftn(opt.p1*compThx4d(v1_k+e11_k,opt)+opt.p1*(v2_k+e12_k)+opt.p3*(opt.M.*x_k-e3_k))./Tinv); 
% elseif opt.lambda1~=0 && opt.lambda3==0
% z_k=ifftn(fftn(opt.p1*compThx4d(v1_k+e11_k,opt)+opt.p3*(opt.M.*x_k-e3_k))./Tinv); 
% elseif opt.lambda3~=0 && opt.lambda1==0
% z_k=(opt.p1*(v2_k+e12_k)+opt.p3*(opt.M.*x_k-e3_k))/(opt.p1+opt.p3); 
% end
%--------------------------------------------------------------------------
%% calculate z(k+1)=[p1*R'R+p3*I]-1*[p1*R'(v+e1)+p3*(x-e3)];
zftemp=opt.p3*(opt.M.*x_k-e3_k);

if opt.lambda1~=0||opt.lambda2~=0
    zftemp=zftemp+opt.p1*compThx4d(v1_k+e11_k,opt);
end

if opt.lambda3~=0
    zftemp=zftemp+opt.p1*(v2_k+e12_k);
end

if opt.xf~=0
    zftemp=zftemp+opt.p1*opt.TempIFT(v3_k+e13_k);
end

z_k=ifftn(fftn(zftemp)./Tinv); 

%% solve v1-subproblem and v2-subproblem
if opt.lambda1~=0||opt.lambda2~=0
    v1_k=shrink1(compTx4d(z_k,opt)-e11_k,TVlambda/opt.p1); 
end

if opt.lambda3~=0
    v2_k=LLRp_v3(z_k-e12_k, opt.B, opt.lambda3/opt.p1, opt.p, opt.overlap, opt.do_shift, opt.useParrallel);     
end

%% calculate x(k+1)=[p2*S'S+p3*I]-1*[p2S'(u-e2)+p3(z-e3)];

x_k=(opt.p2*opt.M.*compShx4d(u_k+e2_k,sMaps,opt)+opt.p3*(z_k+e3_k))./SIinv;

%% calculate u(k+1)=[F'F+p2I]-1*[F'y+p2*(Sx-e2)]; this is the only equation with the actual data kU.  
kU=zeros(opt.size,opt.class);
kU(opt.U)=k;

u_k=opt.IFT((kU+opt.p2*opt.FT(compSx4d(opt.M.*x_k,sMaps,opt)-e2_k))./UIinv);

  
%% update residue variables
% e1=e1-Rz+v
if opt.lambda1~=0||opt.lambda2~=0
    e11_k=e11_k-compTx4d(z_k,opt)+v1_k;
end

if opt.lambda3~=0
    e12_k=e12_k-z_k+v2_k;
end

% e2=e2-Sx+u
e2_k=e2_k-compSx4d(opt.M.*x_k,sMaps,opt)+u_k;

% e3=e3-x+z
e3_k=e3_k-opt.M.*x_k+z_k;

%% computing cost function

k_es=opt.FT(compSx4d(opt.M.*x_k,sMaps,opt));
e_data=kU(opt.U)-k_es(opt.U);

obj11=compTx4d(x_k,opt);
% obj12 = norm_evaluate(x_k,opt.p,opt.lambda3,opt.B,0);
cost1(iter)=sum(abs(e_data(:)).^2);
cost2(iter)=sum(TVlambda.*abs(obj11));
cost3(iter)=cost1(iter)+cost2(iter);

if opt.plot&&mod(iter-1,opt.update)==0
    opt.RI.Iter=iter;
    opt.RI.ReconTime=toc;
subplot(1,4,1);imagesc(cat(1,abs(x_k(:,:,3,1)),abs(x_k(:,:,4,1)),abs(x_k(:,:,5,1))),[0 0.6*max(abs(x_k(:)))]);colormap('gray'); title(num2str(iter));axis equal;drawnow
subplot(1,4,2);plot(double(log10(cost1(1:iter)))); title('data consistency');drawnow
subplot(1,4,3);plot(double(cost2(1:iter)));title('l1-norm'); drawnow
subplot(1,4,4);plot(double(cost3(1:iter)));title('Total cost');drawnow
end


iter=iter+1;

% if(mod(iter,10)==0)
% opt.p1=opt.p1*10;
% end
toc

end

imgR=reshape(x_k,[kx,ky,kz,nt]);
clear x_k;

end
