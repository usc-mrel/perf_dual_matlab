
%   Detailed explanation goes here
function [Y] = svt(X,lambda,p)
%function [Y, evaln] = svt(X,lambda,opts)
[U,S,V] = svd(X,'econ'); s = diag(S);
%s = sign(s).*max(s-lambda,0);

%evaln = sum(abs(s(:)).^(opts.p)./(opts.p))*opts.lambda1;
thres=(lambda).*(s.^(p-1));
        
        s=(s-thres);
        s = s.*(s>0);
       
Y = U*diag(s)*V';


end


