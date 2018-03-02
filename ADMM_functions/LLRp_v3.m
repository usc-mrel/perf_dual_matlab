function [g] = LLRp_v3(img, B, lambda, p, overlap, DO_SHIFT, useParrallel)
%function [g, evaln] = LLRp_v3(img,opt, B,lambda,DO_SHIFT)
% LLRP
%   [res_c, res_t] = LLRP(img,opt) performs a Locally Low Rank Proejction on
%   multicoil images. Code is modified from Joshua Trzasko's CLEAR algorithm,
%   ISMRM 2012. res_c is LLR performed after unfolding the coil dimeension.
%   res_t is LLR performed after unfolding the time dimension
%
%   inputs:
%   img: Multi-frame images with size Nx Ny Nt
%   B:   patch size
%   lambda: regularization parameter
%   p:      ...
%   overlap: Allows blocks to overlap (1) or not (0)
%   DO_SHIFT: Allow random shifting between iteration
%   useParrallel: Use parrallel computing
%
%   authors: Sajan Lingala, Yi Guo, Terrence Jao
%   date: 03-01-2014
%   ref: Trazasko, Joshua Manduca, Armando.CLEAR: Calibration-Free Parallel
%   Imaging using Locally Low-Rank Encouraging Reconstruction. ISMRM 2012.
%   0517

%--- Parameters ---%


if overlap == 1
    Nshift = 8; %Always shift in 8 directions plus 1 without shift. 
else 
    Nshift = 1;
end 

%--- Peform LLR ---%
[Nx,Ny,Nz,Nt] = size(img);

if DO_SHIFT
    vec = randperm(B+1)-B/2-1; dx = round(vec(1));
    vec = randperm(B+1)-B/2-1; dy = round(vec(1));
    img = circshift(img,[dx dy]);
end

g = resCalc(img);

if DO_SHIFT
    g = circshift(g,[-dx -dy]);
end

    function res = resCalc(img)
        
        res = zeros([Nx, Ny, Nz, Nt]);
      
        %Find linear indices of all Patches projected along time frames
        ind = patchInd_4d(img,B,overlap);
        if ~overlap
            nP = numel(ind)/(Nt*B*B*B); %number of patches created
        else
            nP = numel(ind)/(Nt*B*B*B*Nshift);
        end 
        BLOCK = img(ind);
        %keyboard;
        if ~ overlap
            BLOCK = reshape(BLOCK, [B*B*B, nP, Nt]);
            BLOCK = permute(BLOCK, [1 3 2]);
        else 
            BLOCK = reshape(BLOCK, [B*B*B, nP, Nt, Nshift]);
            BLOCK = permute(BLOCK, [1 3 2 4]);
            BLOCK = reshape(BLOCK, [B*B*B, Nt, nP*Nshift]);
        end

        % Parallel Computing Implementation %
        %evaln = zeros(nP*Nshift,1);
        if useParrallel
            
            parfor ip = 1:nP*Nshift
                %[BLOCK(:,:,p) evaln(p)]= svt(BLOCK(:,:,p), lambda,opt);
                [BLOCK(:,:,ip)]= svt(BLOCK(:,:,ip), lambda, p);
            end
            
        
        %evaln = sum(evaln);
        % Array Function Implementation %
        else
            for ip = 1:nP*Nshift
                patchArray(ip).data = BLOCK(:,:,ip);
            end
            %evaln = 0;
            svtFunc = @shrink;
            % PatchArray is seen by shrink, and being modified within shrink. 
            % No need to give arrayfun an  output
            arrayfun(svtFunc, patchArray, 'UniformOutput', false);
            BLOCK = cat(3,patchArray.data);

        end
        if ~ overlap
            BLOCK = permute(BLOCK, [1 3 2]);
        else 
            BLOCK = reshape(BLOCK, [B*B*B, Nt, nP, Nshift]);
            BLOCK = permute(BLOCK, [1 3 2 4]);
        end
        
        tempres = res;
        for k = 1:Nshift
            tempres(ind(:,k)) = BLOCK(:,:,:,k);
            res = res + tempres;
        end 

        %evaln = evaln/Nshift;
        res = res/Nshift; %could be incorrect for some regions that did not overlap Nshift times...

    end

    function patch = shrink(patch)
    %    %keyboard;
        [U,S,V]=svd(patch.data, 'econ');
        s = diag(S);
    %    %s = sign(s).*max(s-lambda,0);
        %evaln = evaln + sum(abs(s(:)).^(opt.p)./(opt.p))*opt.lambda1;
        thres=(lambda).*(s.^(p-1));
    %    %thres=(lambda/beta).*(s.^(opt.p-1));
        s=(s-thres);
        s = s.*(s>0);
        patch = U*diag(s)*V';
        
    end
end
