function Thx=compThx3d(x,opt)
% Calculate 2d spatial transposed Total variation for 3d data set

begin=1;

kx=opt.size(1);
ky=opt.size(2);
nt=opt.size(3);
NL=kx*ky*nt;

TV1=reshape(x(begin:begin+NL-1),kx,ky,nt); % split the vecotr back to TV coefficient
begin=begin+NL;

TV2=reshape(x(begin:begin+NL-1),kx,ky,nt);
begin=begin+NL;

TV3=reshape(x(begin:begin+NL-1),kx,ky,nt);
begin=begin+NL;

TVh1=TV1-circshift(TV1,[-1 0 0]);
TVh2=TV2-circshift(TV2,[0 -1 0]);
TVh3=TV3-circshift(TV3,[0 0 -1]);

Thx=TVh1+TVh2+TVh3;
end