function imout = imresize_perf(imin,usercv,diastole)

    if diastole
        opxres = usercv(7);
        opyres = usercv(8);
        opzres = usercv(9);
    else
        opxres = usercv(10);
        opyres = usercv(11);
        opzres = usercv(12);
    end
    
    imout = zeros(opxres,opxres,opzres,size(imin,4));
      
    for m=1:size(imin,4)
        for n=1:opzres
            imout(:,:,n,m) = imresize(imin(:,:,n,m),[opxres, opxres]);
        end
    end
    

end