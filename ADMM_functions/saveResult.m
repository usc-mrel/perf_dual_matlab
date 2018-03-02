function saveResult(X_est,subdirectory,opts)
%---------------------------------------------------------------
% Write the output
%---------------------------------------------------------------
%Create Results Directory if Missing 
system(sprintf('mkdir -p %s', subdirectory));
cd(sprintf('%s', subdirectory));

%Write the output

fn1 = sprintf('lam1_%0.2f_lam2_%0.2f_lam3_%0.2f_patch_%d_Olap_%d_.mat',opts.lambda1, opts.lambda2, opts.lambda3, opts.B, opts.overlap);
save(fn1,'X_est','opts.nerror');
cd(sprintf('../..'));

end
