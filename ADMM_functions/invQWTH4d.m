function Tinv=invQWTH4d(opt)
% function to calculate inverse of [R'R+p*I] directly
% R=[W;T], is 4d wavelet,3d TV and temporal filter
% Yi Guo 05/17/2013

% to 4d 05/20

kx=opt.size(1);
ky=opt.size(2);
kz=opt.size(3);
nt=opt.size(4);

if opt.beta1~=0
    p1=opt.p1;
else
    p1=0;
end


t1=[1;-1];

if opt.beta2~=0
   
T1z=zeros(kx,ky,kz,nt);
T2z=T1z;
T3z=T1z;

% three direction TV
T1z(1:2,1,1,1)=t1;
T2z(1,1:2,1,1)=t1;
T3z(1,1,1:2,1)=t1;

T1=fftn(T1z);
T2=fftn(T2z);
T3=fftn(T3z);
if opt.size(3)==1
    T3=0;
end

else 
    T1=0;
    T2=0;
    T3=0;
end


% make phase shift in Fourier domain to move circushift inside
[ph2,ph1,ph3]=meshgrid(0:ky-1,0:kx-1,0:kz-1);

ph1=repmat(ph1,[1 1 1 nt]);
ph2=repmat(ph2,[1 1 1 nt]);
ph3=repmat(ph3,[1 1 1 nt]);

[~,~,ph4]=meshgrid(0:kz-1,0:ky-1,0:nt-1);
ph4=repmat(ph4,[1 1 1 kx]);
ph4=permute(ph4,[4 1 2 3]);

shift1=exp(1j*2*pi*ph1/kx);
shift2=exp(1j*2*pi*ph2/ky);
shift3=exp(1j*2*pi*ph3/kz);
shift4=exp(1j*2*pi*ph4/nt);

if opt.beta3~=0&&opt.size(4)>5
HH=zeros(kx,ky,kz,nt);
% HH(1,1,1,1:5)=[0.06 0.44 -1 0.44 0.06];
HH(1,1,1,1:2)=[1;-1];
HHF=fftn(HH);
%W=cat(4,zeros(kx,ky,kz,1),ones(kx,ky,kz,5)*1,ones(kx,ky,kz,3)*0.1,ones(kx,ky,kz,26)*2);

else 
    HHF=0;
end

    Tinv=opt.p1*(-T1.*shift1.*T1-T2.*shift2.*T2-T3.*shift3.*T3-HHF.*shift4.*HHF)+p1+opt.p3;
end

