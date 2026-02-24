function result = UBt(I,hsize)

    [rows, cols, chan] = size(I);

    % ガウスカーネルによるブラーリングの転置
    psf = fspecial('average',hsize);
    psfsize = size(psf);
    blu = zeros(rows, cols);
    blu(1:psfsize(1), 1:psfsize(2)) = psf;
    blu = gpuArray(circshift(blu, [-ceil((psfsize-1)/2) -ceil((psfsize-1)/2)]));
    bluf = fft2(blu);
    bluft = conj(bluf);
    bluft = repmat(bluft, [1,1,chan]);
    result = real(ifft2((fft2(I)).*bluft));
end