% input : Hr, Lr, Lt

%% Initial values
preHr  = gpuArray(ones(size(Hr))); % Estimate of Noiseless Reference HR
preHt  = gpuArray(ones(size(Hr))); % Estimate of Noiseless Target HR
Z1 = gpuArray(ones(size(WD(Hr))));
Z2 = gpuArray(ones(size(WD(Hr))));
Z3 = gpuArray(ones(size(WD(Hr))));
Z4 = gpuArray(ones(size(Hr)));
Z5 = gpuArray(ones(size(Lr)));
Z6 = gpuArray(ones(size(Lt)));

objfuncval_list  = zeros(1, params.max_iteration);
psnr_Hr_list     = zeros(1,params.max_iteration);
psnr_Ht_list     = zeros(1,params.max_iteration);
error_preHr_list = zeros(1,params.max_iteration);
error_preHt_list = zeros(1,params.max_iteration);
alpha_list       = zeros(1,params.max_iteration);
error_alpha_list = zeros(1,params.max_iteration);
r1_list          = zeros(1, params.max_iteration);
r2_list          = zeros(1, params.max_iteration);
mean_dif_list    = zeros(1, params.max_iteration);
ite_time_list    = zeros(1, params.max_iteration);

stlm_r = stretchlim(Hr_GT(:, :,[4 3 2]));
stlm_t = stretchlim(Ht_GT(:, :,[4 3 2]));
stlm = [min([stlm_r; stlm_t]); max([stlm_r; stlm_t])];

%% Optimization
is_converged = 0;
iteration = 1;
while is_converged == 0
    
    tic;

    % Update Primal Variables
    Ur = DtWt(Z1) + DtWt(Z3) + Z4 + BtSt(Z5);   % line 4
    Ut = DtWt(Z2) - DtWt(Z3) + BtSt(Z6);        % line 5
    preHr_next = preHr - params.gamma1(1) * Ur; % line 6
    preHt_next = preHt - params.gamma1(2) * Ut; % line 7
    for chan = 1:chans
        preHr_next(:,:,chan) = indicator_hyperslab(preHr_next(:,:,chan), ...
                                                   params.lowr(1,1,chan), ...
                                                   params.highr(1,1,chan)); % line 8
        preHt_next(:,:,chan) = indicator_hyperslab(preHt_next(:,:,chan), ...
                                                   params.lowt(1,1,chan), ...
                                                   params.hight(1,1,chan)); % line 9
    end

    % Over-relaxed variables -- lines 16-19
    Hr_bar  = 2*preHr_next - preHr; 
    Ht_bar  = 2*preHt_next - preHt; 

    % Update parameter alpha -- line 20
    params.alpha_next = params.alpha_coef * LqNorm(WD(preHr_next)) * Dif_L;
    alpha_list(iteration) = params.alpha_next;
    
    % Update dual variables
    Z1_next = Z1 + params.gamma2(1)*WD(Hr_bar); % line 21
    Z2_next = Z2 + params.gamma2(2)*WD(Ht_bar); % line 22
    Z3_next = Z3 + params.gamma2(3)*(WD(Hr_bar)-WD(Ht_bar)); % line 23
    Z4_next = Z4 + params.gamma2(4)*Hr_bar; % line 24
    Z5_next = Z5 + params.gamma2(5)*SB(Hr_bar); % line 25
    Z6_next = Z6 + params.gamma2(6)*SB(Ht_bar); % line 26
    
    Z1_next = Z1_next - params.gamma2(1)*prox12band(Z1_next/params.gamma2(1), 1/params.gamma2(1)); % line 30
    Z2_next = Z2_next - params.gamma2(2)*prox12band(Z2_next/params.gamma2(2), params.lambda/params.gamma2(2)); % line 31
    Z3_next = Z3_next - params.gamma2(3)*projLqBall(Z3_next/params.gamma2(3), params.alpha_next); % line 32
    Z4_next = Z4_next - params.gamma2(4)*proj_Fball(Z4_next/params.gamma2(4), Hr, params.epsilonh); % line 33
    Z5_next = Z5_next - params.gamma2(5)*proj_Fball(Z5_next/params.gamma2(5), Lr, params.epsilonl); % line 34
    Z6_next = Z6_next - params.gamma2(6)*proj_Fball(Z6_next/params.gamma2(6), Lt, params.epsilonl); % line 35

    ite_time = toc;
    ite_time_list(iteration) = ite_time;

    % Compute Error
    error_preHr = (sqrt(sum((preHr - preHr_next).^2, "all")/sum(preHr.^2, "all")));
    error_preHt = (sqrt(sum((preHt - preHt_next).^2, "all")/sum(preHt.^2, "all")));
    error_alpha = abs(params.alpha - params.alpha_next)/abs(params.alpha);
    error_preHr_list(iteration) = error_preHr;
    error_preHt_list(iteration) = error_preHt;
    error_alpha_list(iteration) = error_alpha;
    
    r1 = sqrt(sum((Lr - SB(preHr)).^2,'all'));
    r2 = sqrt(sum((Lt - SB(preHt)).^2,'all'));

    % Check Convergence
    if ((error_preHr <= params.stopping_criterion && ...
         error_preHt <= params.stopping_criterion && ...
         r1 <= params.epsilonl && ...
         r2 <= params.epsilonl) || ...
         iteration >= params.max_iteration)
        is_converged = 1;
    end
    
    preHr = preHr_next;
    preHt = preHt_next;
    params.alpha = params.alpha_next;
    Z1 = Z1_next;
    Z2 = Z2_next;
    Z3 = Z3_next;
    Z4 = Z4_next;
    Z5 = Z5_next;
    Z6 = Z6_next;
    
    objfuncval = l12norm(WD(preHr)) + params.lambda*l12norm(WD(preHt));
    [psnr_preHr,~,~,~,~]  = cal_metrics(gather(preHr),Hr_GT);
    [psnr_preHt,~,~,~,~]  = cal_metrics(gather(preHt),Ht_GT);
    rmse = sqrt(immse(gather(preHt), Ht_GT));
    mean_dif = mean(Lt(:,:,1),'all') - mean(preHt(:,:,1),'all');

    objfuncval_list(iteration) = objfuncval;
    psnr_Hr_list(iteration) = psnr_preHr;
    psnr_Ht_list(iteration) = psnr_preHt;
    r1_list(iteration) = r1;
    r2_list(iteration) = r2;
    mean_dif_list(iteration) = mean_dif;

    if mod(iteration,100) == 0
        
        disp("iter. : " + num2str(iteration,'%05d') + ...
            ", objval : " + num2str(objfuncval,'%.2f') + ...
            ", PSNR of Hr : " + num2str(psnr_preHr,'%.2f') + ...
            ", PSNR of Ht : " + num2str(psnr_preHt,'%.2f') + ...
            ", RMSE of Ht : " + num2str(rmse,'%.4f') + ...
            ", iter. time : "+ num2str(sum(ite_time_list)))
        
        figure(1);
        SBHr = SB(preHr);
        SBHt = SB(preHt);
        
        subplot(2,4,1), imshow(imadjust(Hr(:,:,band_for_show),stlm,[])), title('$$ \mathbf{h_r} $$','Interpreter','latex','FontSize',20)
        subplot(2,4,2), imshow(imadjust(Lr(:,:,band_for_show),stlm,[])), title('$$ \mathbf{l_r} $$','Interpreter','latex','FontSize',20)
        subplot(2,4,3), imshow(imadjust(Ht_GT(:,:,band_for_show),stlm,[])), title('$$ (\mathbf{\hat{h}_t}) $$','Interpreter','latex','FontSize',20)
        subplot(2,4,4), imshow(imadjust(Lt(:,:,band_for_show),stlm,[])), title('$$ \mathbf{l_t} $$','Interpreter','latex','FontSize',20)
        subplot(2,4,5), imshow(imadjust(preHr(:,:,band_for_show),stlm,[])), title('$$ \mathbf{\tilde{h}_r} $$','Interpreter','latex','FontSize',20)
        subplot(2,4,6), imshow(imadjust(SBHr(:,:,band_for_show),stlm,[])), title('$$ \mathbf{SB\tilde{h}_r} $$','Interpreter','latex','FontSize',20)
        subplot(2,4,7), imshow(imadjust(preHt(:,:,band_for_show),stlm,[])), title('$$ \mathbf{\tilde{h}_t} $$','Interpreter','latex','FontSize',20)
        subplot(2,4,8), imshow(imadjust(SBHt(:,:,band_for_show),stlm,[])), title('$$ \mathbf{SB\tilde{h}_t} $$','Interpreter','latex','FontSize',20)
        saveas(gcf,append(output_dir,'/current_result.png'))

        figure(2);
        subplot(1,2,1), imshow(imadjust(preHt(90:160,240:310,band_for_show),stlm,[]))
        subplot(1,2,2), imshow(imadjust(Ht_GT(90:160,240:310,band_for_show),stlm,[]))

        figure(3);
        subplot(3,1,1), plot(objfuncval_list(1:iteration)), title('Objective Function Value');
        subplot(3,1,2), plot(psnr_Ht_list(1:iteration)), title('PSNR of Ht');
        subplot(3,1,3), plot(alpha_list(1:iteration)), title('Alpha');
    end
    iteration = iteration + 1;
end

preHr = gather(preHr);
preHt = gather(preHt);

save(output_file,...
    'preHr','preHt',...
    'objfuncval_list','psnr_Hr_list','psnr_Ht_list',...
    'mean_dif_list','ite_time_list',...
    'error_preHr_list','error_preHt_list', 'error_alpha_list',...
    'r1_list','r2_list', 'alpha_list')
save(alpha_file, 'alpha_list')
save(preHt_file, 'preHt')