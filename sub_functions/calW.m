function [Wv, Wh, Wbr, Wbl] = calW(Hr,k,delta,pre_denoise)
    
    if pre_denoise == 1
        % medianフィルタで pre noise removal
        med_Hr = Hr;
        for b = 1:size(Hr,3)
            med_Hr(:,:,b) = medfilt2(Hr(:,:,b),[3 3]);
        end
        Hr_band_mean = mean(med_Hr,3);
    else
        Hr_band_mean = mean(Hr,3);
    end
    
    [D, Dt] = GenerateDifferenceOperator();
    DHr_band_mean = D(Hr_band_mean);

    Wv = exp(-(DHr_band_mean(:,:,:,1)./delta).^2);
    Wh = exp(-(DHr_band_mean(:,:,:,2)./delta).^2);
    Wbr = exp(-(DHr_band_mean(:,:,:,3)./delta).^2);
    Wbl = exp(-(DHr_band_mean(:,:,:,4)./delta).^2);

    W_mat = cat(3, Wv, Wh, Wbr, Wbl);
    W_mat = filter_top_k(W_mat, k);

    Wv = W_mat(:,:,1);
    Wh = W_mat(:,:,2);
    Wbr = W_mat(:,:,3);
    Wbl = W_mat(:,:,4);
end

