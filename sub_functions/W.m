function [result] = W(Dz, Wv, Wh, Wbr, Wbl)
    
    chans = size(Dz,3);
    Wv = repmat(Wv, [1,1,chans]);

    Dvz = Dz(:,:,:,1);
    Dhz = Dz(:,:,:,2);
    Dbrz = Dz(:,:,:,3);
    Dblz = Dz(:,:,:,4);
    
    WvDvz = Wv.*Dvz;
    WhDhz = Wh.*Dhz;
    WbrDbrz = Wbr.*Dbrz;
    WblDblz = Wbl.*Dblz;

    result = cat(4, WvDvz, WhDhz, WbrDbrz, WblDblz);
end

