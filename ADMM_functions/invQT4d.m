function Tinv=invQT4d(opt)

if opt.lambda3~=0
    p1=opt.p1;
else
    p1=0;
end

kx=opt.size(1);
ky=opt.size(2);
kz=opt.size(3);
nt=opt.size(4);

T1z=zeros(kx,ky,kz,nt);
T2z=T1z;
T3z=T2z;
T4z=T3z;

t1=[1;-1];

if (opt.lambda1==0)
    T1=0;
    T2=0;
    T3=0;
else
    % three direction TV
    T1z(1:2,1,1,1)=t1;
    T2z(1,1:2,1,1)=t1;
    T3z(1,1,1:2,1)=t1;

    T1=fftn(T1z);
    T2=fftn(T2z);
    T3=fftn(T3z);

end

if (opt.lambda2==0)
    T4=0;
else

    T4z(1,1,1,1:2)=t1;

    T4=fftn(T4z);
end


% make phase shift in Fourier domain to move circushift inside
[ph1,ph2,ph3,ph4]=ndgrid(0:kx-1,0:ky-1,0:kz-1,0:nt-1);

shift1=exp(1j*2*pi*ph1/kx);
shift2=exp(1j*2*pi*ph2/ky);
shift3=exp(1j*2*pi*ph3/kz);
shift4=exp(1j*2*pi*ph4/nt);

Tinv=opt.p1*(-T1.*shift1.*T1-T2.*shift2.*T2-T3.*shift3.*T3-T4.*shift4.*T4)+p1+opt.p3;

