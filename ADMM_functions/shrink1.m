function xres=shrink1(x,lambda)
% shrink function to approximate l1 norm
% Yi Guo 04/29/2013

xres=sign(x).*max(abs(x)-lambda,0);
end