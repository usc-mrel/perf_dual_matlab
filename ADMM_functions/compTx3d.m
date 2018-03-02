function [Tx,LL]=compTx3d(x,opt)
% Calculate 2d spatial Total variation for 3d data set

TV1=x-circshift(x,[1 0 0]);
TV2=x-circshift(x,[0 1 0]);
TV3=x-circshift(x,[0 0 1]);
Tx=[TV1(:);TV2(:);TV3(:)];
LL=numel(x)*3;
end