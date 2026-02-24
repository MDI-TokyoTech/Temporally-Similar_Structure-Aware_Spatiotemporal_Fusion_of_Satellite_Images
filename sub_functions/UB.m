function result = UB(I,hsize)
    [rows, cols, chan] = size(I);

    psf = fspecial('average',hsize);
    psfsize = size(psf);
    blu = zeros(rows, cols);
    blu(1:psfsize(1), 1:psfsize(2)) = psf;
    blu = gpuArray(circshift(blu, [-ceil((psfsize-1)/2) -ceil((psfsize-1)/2)]));
    bluf = fft2(blu);
    bluf = repmat(bluf, [1,1,chan]);
    result = real(ifft2((fft2(I)).*bluf));    
end