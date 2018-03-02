function [kd_zd, idxd, ks_zd, idxs,usercv] = sampling(dr,rhfrsize,usercv,rhnframes,rhnslices,ncoil, accel_ky, accel_kz)

nsegs = usercv(3);
nsegs_accel = usercv(4);

opxres = usercv(7);
opyres = usercv(8);
opzres = usercv(9);

opxres_sys = usercv(10);
opyres_sys = usercv(11);
opzres_sys = usercv(12);

skedtype = usercv(14);
skedtype_sys = usercv(15);
sked_idx = usercv(16);
sked_idx_sys = usercv(17);

nsegs_sys = usercv(18);
nsegs_accel_sys = usercv(19);

if skedtype == 2
    switch sked_idx
        case 0 
            skedfile = 'SKED/skedR9_mtx116_76_12_fov280_120_r4'; 
            opzres = 12; %hack!
            usercv(9) = opzres;
        case 1
            skedfile = 'SKED/skedR9_mtx78_52_10_fov280_100_r3';
        case 2
            skedfile = 'SKED/skedR9_mtx116_76_10_fov280_100_r4';
            opzres = 10; %hack!
            usercv(9) = opzres;
        case 3
            skedfile = 'SKED/skedR9_mtx78_52_8_fov280_80_r3';
    end
end

if skedtype_sys == 2
    switch sked_idx_sys
        case 0 
            skedfile_sys = 'SKED/skedR9_mtx116_76_12_fov280_120_r4'; 
        case 1
            skedfile_sys = 'SKED/skedR9_mtx78_52_10_fov280_100_r3';
            opzres_sys = 10; %hack!
            usercv(12) = opzres_sys;
        case 2
            skedfile_sys = 'SKED/skedR9_mtx116_76_10_fov280_100_r4';
        case 3
            skedfile_sys = 'SKED/skedR9_mtx78_52_8_fov280_80_r3';
            opzres_sys = 8; %hack!
            usercv(12) = opzres_sys;
    end
    
end

%% sort k data in diastole and systole
% Dimensions:
% 1: rhfrsize = MAX( opxres, opxres_sys );
% 2: rhnframes+1 = 1+ MAX(nsegs_accel,nsegs_accel_sys);
% 3: rhnecho = 2
% 4: rhnslices = nframe;
% 5: ncoils


dr = reshape(dr,[rhfrsize, (rhnframes+1), 2, rhnslices, ncoil]);
kd = squeeze(dr(:, 2:(rhnframes+1), 1, :, :)); 
ks = squeeze(dr(:, 2:(rhnframes+1), 2, :, :));

if (rhfrsize > opxres) %rhnfrsize = opxres_sys > opxres
    kd = kd(1:opxres,:,:,:);
else %rhnfrsize = opxres > opxres_sys
    ks = ks(1:opxres_sys,:,:,:);
end

if (rhnframes > nsegs_accel)
    kd = kd(:,1:nsegs_accel,:,:);
else
    ks = ks(:,1:nsegs_accel_sys,:,:);
end

nframe = rhnslices;

%% zeros padding recon
% sampling pattern repeat itself from 1--> nsegs/nsegs_accel

kd_zd = zeros(opxres, opyres, opzres, nframe, ncoil);
idxd = zeros(opyres,opzres,nframe);

%---for kd_zd
if (skedtype==1)
    disp('regular undersample');
    %-----
    [ mask sampidx] = mask_lattice(opxres, opyres, opzres,accel_ky,accel_kz );
    %-----
    
    sked_ky = real(sampidx(:));
    sked_kz = imag(sampidx(:));
    
    for nnt=1:nframe
        for nnp = 1:nsegs_accel
            
            tmp = mod(nnt,accel_ky*accel_kz);
            if tmp ==0
                tmp = accel_ky*accel_kz;
            end
            
            nn = (tmp-1)*nsegs_accel+nnp;
            
            kd_zd(:,sked_ky(nn), sked_kz(nn), nnt,:) = reshape(kd(:,nnp,nnt,:),[opxres 1 1 1 ncoil]);
            idxd(sked_ky(nn), sked_kz(nn),nnt) = idxd(sked_ky(nn), sked_kz(nn),nnt)+1;
        end
        
    end
    
elseif (skedtype==2)
    disp('irregular undersample');
    disp(['Sked file index:' num2str(sked_idx)]);
    
    [sked_ky,sked_kz] = readSkedFile(skedfile);
    
    for nnt=1:nframe
        for nnp = 1:nsegs_accel
            
            tmp = mod(nnt,nsegs/nsegs_accel);
            if tmp ==0
                tmp = nsegs/nsegs_accel;
            end
            
            nn = (tmp-1)*nsegs_accel+nnp;
            
            kd_zd(:,sked_ky(nn)+1, sked_kz(nn)+1, nnt,:) = reshape(kd(:,nnp,nnt,:),[opxres 1 1 1 ncoil]);
            
            disp(['nnt=' num2str(nnt) 'nnp=' num2str(nnp)]);
            idxd(sked_ky(nn)+1, sked_kz(nn)+1,nnt) = idxd(sked_ky(nn)+1, sked_kz(nn)+1,nnt)+1;

        end
    end
    
end

%%

ks_zd = zeros(opxres_sys, opyres_sys, opzres_sys, nframe, ncoil);
idxs = zeros(opyres_sys,opzres_sys,nframe);

%---for ks_zd
if (skedtype_sys==1)
    disp('regular undersample');
    
    [ mask_sys sampidx_sys] = mask_lattice(opxres_sys, opyres_sys, opzres_sys, accel_ky,accel_kz );
    
    sked_ky = real(sampidx_sys(:));
    sked_kz = imag(sampidx_sys(:));
    
    for nnt=1:nframe
        for nnp = 1:nsegs_accel_sys
            
            tmp = mod(nnt,accel_ky*accel_kz);
            if tmp ==0
                tmp = accel_ky*accel_kz;
            end
            
            nn = (tmp-1)*nsegs_accel_sys+nnp;
            
            ks_zd(:,sked_ky(nn), sked_kz(nn),nnt, :) = reshape(ks(:,nnp,nnt,:),[opxres_sys 1 1 1 ncoil]);
            idxs(sked_ky(nn), sked_kz(nn),nnt) = idxs(sked_ky(nn), sked_kz(nn),nnt)+1;
        end
    end
    
elseif (skedtype_sys==2)
    disp('irregular undersample');
    disp(['Sked file index:' num2str(sked_idx_sys)]);
    
    [sked_ky,sked_kz] = readSkedFile(skedfile_sys);
    
    for nnt=1:nframe
        for nnp = 1:nsegs_accel_sys
            
            tmp = mod(nnt,nsegs_sys/nsegs_accel_sys);
            if tmp ==0
                tmp =nsegs_sys/nsegs_accel_sys;
            end
            
            nn = (tmp-1)*nsegs_accel_sys+nnp;
            
            ks_zd(:,sked_ky(nn)+1, sked_kz(nn)+1, nnt,:) = reshape(ks(:,nnp,nnt,:),[opxres_sys 1 1 1 ncoil]);
            idxs(sked_ky(nn)+1, sked_kz(nn)+1,nnt) = idxs(sked_ky(nn)+1, sked_kz(nn)+1,nnt)+1;

        end
    end
    
end

end