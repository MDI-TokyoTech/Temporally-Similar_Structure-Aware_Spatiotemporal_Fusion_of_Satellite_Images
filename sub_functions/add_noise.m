function [Hr, Ht, Lr, Lt] = add_noise(Hr_GT, Ht_GT, Lr_GT, Lt_GT, params)

    Hr_GT(Hr_GT<0) = 0; Hr_GT(Hr_GT>1) = 1;
    Lr_GT(Lr_GT<0) = 0; Lr_GT(Lr_GT>1) = 1;
    Ht_GT(Ht_GT<0) = 0; Ht_GT(Ht_GT>1) = 1;
    Lt_GT(Lt_GT<0) = 0; Lt_GT(Lt_GT>1) = 1;

    deg_HR = struct();
    deg_HR.scaling_factor = params.HR_poisson;
    deg_HR.stripe_rate = params.HR_stripe_rate;
    deg_HR.stripe_intensity = params.HR_stripe_intensity;
    deg_HR.Gaussian_sigma = params.HR_gaussian;
    deg_HR.sparse_rate = params.HR_sparse;

    deg_LR = struct();
    deg_LR.scaling_factor = params.LR_poisson;
    deg_LR.stripe_rate = params.LR_stripe_rate;
    deg_LR.stripe_intensity = params.LR_stripe_intensity;
    deg_LR.Gaussian_sigma = params.LR_gaussian;
    deg_LR.sparse_rate = params.LR_sparse;

    Hr = Generate_obsv(Hr_GT, deg_HR);
    Ht = Generate_obsv(Ht_GT, deg_HR);
    Lr = Generate_obsv(Lr_GT, deg_LR);
    Lt = Generate_obsv(Lt_GT, deg_LR);

    Hr(Hr<0) = 0; Hr(Hr>1) = 1;
    Ht(Ht<0) = 0; Ht(Ht>1) = 1;
    Lr(Lr<0) = 0; Lr(Lr>1) = 1;
    Lt(Lt<0) = 0; Lt(Lt>1) = 1;

end

function I_noisy = Generate_obsv(I_clean, deg)
    [n1, n2, n3] = size(I_clean);
    %% Generating poisson noise
    scaling_factor = deg.scaling_factor;
    if scaling_factor == 0
        I_noisy = I_clean;
    else
        I_noisy = poissrnd(I_clean*scaling_factor)/scaling_factor;
    end

    %% Generating stripe noise
    stripe_rate             = deg.stripe_rate;
    if stripe_rate > 0
        stripe_intensity        = deg.stripe_intensity;
        
        sparse_stripe = 2*(imnoise(0.5*ones(1, n2, n3), "salt & pepper", stripe_rate) - 0.5).* ...
            rand(1, n2, n3).*ones(n1, n2, n3);
        if max(abs(sparse_stripe), [], "all") ~= 0
            stripe_noise = stripe_intensity.*sparse_stripe./max(abs(sparse_stripe), [], "all");
            I_noisy = I_noisy + stripe_noise;
        else
            stripe_noise = zeros(n1, n2, n3);
        end
        
        deg.stripe_noise = stripe_noise;
    end
    
    
    %% Generating Gaussian noise
    Gaussian_sigma = deg.Gaussian_sigma;
    
    Gaussian_noise = Gaussian_sigma*randn(n1, n2, n3);
    
    I_noisy = I_noisy + Gaussian_noise;
    deg.Gaussian_noise = Gaussian_noise;
    
    
    %% Generating sparse noise
    sparse_rate = deg.sparse_rate;
    
    HSI_tmp = I_noisy;
    
    Sp = 0.5*ones(n1, n2, n3);
    Sp = imnoise(Sp,'salt & pepper',sparse_rate);
    
    I_noisy(Sp==0) = 0;
    I_noisy(Sp==1) = 1;
    
    deg.sparse_noise = I_noisy - HSI_tmp;
end