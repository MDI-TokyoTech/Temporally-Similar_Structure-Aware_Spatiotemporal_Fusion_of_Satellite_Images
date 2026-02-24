addpath(genpath('./sub_functions'));

%% HyperParameter Setting
params.q = 12;         % Norm q in TGEC. You can choose 12(the mixed l12-norm), 2(the l2-norm), or 1(the l1-norm). Recommended 12.
params.alpha_coef = 5; % Coefficent for TGEC. Recommended 5.
params.delta = 0.1;    % Parameter to control the distribution of weights. Recommended 0.1.
params.k = 2;          % Parameter to specify how many of the four directional weights are forced to be zero. Recommended 2.
params.lambda = 1;     % Balance parameter in (9).
params.spatial_resolution_gap = 20; % Spatial resolution gap between HR and LR.
params.max_iteration = 10000; % Maximum number of algorithm iterations.

%% Noise Setting
params.HR_gaussian         = 0; % Standard Deviation of Gaussian noise in HR image.
params.HR_sparse           = 0; % Superimposition rate of sparse noise in HR image.
params.HR_stripe_rate      = 0.05; % Superimposition rate of stripe noise in HR image.
params.HR_stripe_intensity = 0.2; % Intensity of stripe noise in HR image.
params.HR_poisson          = 0; % Scaling factor of Poisson noise in HR image.
params.LR_gaussian         = 0; % Standard Deviation of Gaussian noise in LR image.
params.LR_sparse           = 0; % Superimposition rate of sparse noise in LR image.
params.LR_stripe_rate      = 0.05; % Superimposition rate of stripe noise in LR image.
params.LR_stripe_intensity = 0.2; % Intensity of stripe noise in LR image.
params.LR_poisson          = 0; % Scaling factor of Poisson noise in LR image.
params.stopping_criterion  = 10^(-5); % The algotithm stops when the update error is less than this value.
% Pre-denoising is needed or not.
if params.HR_poisson == 0 && params.HR_gaussian == 0 && params.HR_sparse == 0 && params.HR_stripe_rate == 0
    params.pre_denoise = 0;
else
    params.pre_denoise = 1;
end

%% Load Dataset and Add Noise
image_file = 'dataset/sample_data.mat';
load(image_file) % load Hr_GT, Ht_GT, Lr_GT, Lt_GT
nH = numel(Hr_GT);
nL = numel(Lr_GT);
[rowsH, colsH, chans] = size(Hr_GT);
[rowsL, colsL, chans] = size(Lr_GT);
nHc = rowsH*colsH; 
nLc = rowsL*colsL;
[Hr, Ht, Lr, Lt] = add_noise(Hr_GT, Ht_GT, Lr_GT, Lt_GT, params);

%% 
output_dir = './Results/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_file  = append(output_dir,'/final_result.mat');
alpha_file   = append(output_dir,'/alpha_transitions.mat');
preHt_file   = append(output_dir,'/preHt_result.mat');
obsdata_file = append(output_dir,'/obsdata.mat');
save(obsdata_file,'Hr','Ht','Lr','Lt');


band_for_show = [4 3 2];

%% Define Operators
% Blurring and Downsampling Operators
B = @(z) UB(z,params.spatial_resolution_gap);
Bt = @(z) UBt(z,params.spatial_resolution_gap);
SB = @(z) S(B(z),params.spatial_resolution_gap,'c');
BtSt = @(z) Bt(St(z,params.spatial_resolution_gap,'c'));
% Difference Operator and Weight Matrix
[D, Dt] = GenerateDifferenceOperator();
[Wv, Wh, Wbr, Wbl] = calW(Hr,params.k,params.delta,params.pre_denoise);
WD = @(z) W(D(z), Wv, Wh, Wbr, Wbl);
DtWt = @(z) Dt(W(z, Wv, Wh, Wbr, Wbl));
W_op = max(cat(3,Wv,Wh,Wbr,Wbl),[],'all'); % Operator Norm of Weight Matrix
% Lq-norm and Projection onto Lp Ball
if params.q == 1
    projLqBall = @(z,alpha) projFastL1Ball(z, alpha);
    LqNorm = @(z) sum(abs(z),"all");
elseif params.q == 2
    projLqBall = @(z,alpha) proj_Fball(z, zeros(size(z)), alpha);
    LqNorm = @(z) sqrt(sum(z.^2,"all"));
elseif params.q == 12
    projLqBall = @(z,alpha) proj_L12(z, alpha, [3,4]);
    LqNorm = @(z) sum(sqrt(sum(z.^2,[3,4])),'all'); % z = W×H×B×4
end

%% Formulation Parameter Setting
% Data-Fidelity Constraint
if params.HR_poisson > 0 % When Poisson noise is added
    params.epsilonh = 0.98 * sqrt((mean(Hr,'all')/params.HR_poisson + params.HR_gaussian^2)*nH*(1-params.HR_sparse));
else % When Poisson noise is not added
    params.epsilonh = 0.98 * params.HR_gaussian * sqrt(nH*(1-params.HR_sparse));
end
params.epsilonl = sqrt(sum((Lr - SB(Hr)).^2,'all'));  % || Lr - SBHr ||_2
params.epsilonlr = params.epsilonl;
params.epsilonlt = params.epsilonl;
% Sparse Noise Constraint
params.etah = 0.98 * params.HR_sparse * nH * 0.5;
params.etal = 0.98 * params.LR_sparse * nL * 0.5;
% Stripe Noise Constraint
params.thetah = 0.98 * params.HR_stripe_intensity * nH * params.HR_stripe_rate * (1- params.HR_sparse)/2;
params.thetal = 0.98 * params.LR_stripe_intensity * nL * params.LR_stripe_rate * (1- params.LR_sparse)/2;

% Edge Constraint
Dif_L = Lr - Lt; % Temporal change
Dif_L = sum(abs(Dif_L(:)))/numel(Lr_GT);
params.alpha = params.alpha_coef * LqNorm(WD(Hr)) * Dif_L;

% brightness constraint
params.beta = zeros(1,1,chans);
params.lowr = zeros(1,1,chans);
params.highr = zeros(1,1,chans);
params.lowt = zeros(1,1,chans);
params.hight = zeros(1,1,chans);
for chan = 1 : chans
    mLr = mean(Lr(:,:,chan),'all');
    mHr = mean(Hr(:,:,chan),'all');
    mLt = mean(Lt(:,:,chan),'all');
    params.beta(1,1,chan) = abs(mLr - mHr); %% 追加
    params.lowr(1,1,chan) = (mLr - params.beta(1,1,chan))*nHc;
    params.highr(1,1,chan) = (mLr + params.beta(1,1,chan))*nHc;
    params.lowt(1,1,chan) = (mLt - params.beta(1,1,chan))*nHc;
    params.hight(1,1,chan) = (mLt + params.beta(1,1,chan))*nHc;
end

save(append(output_dir,'/config.mat'),'-struct','params')
params.gamma1 = [1/(2+32*W_op^2), 1/(1+32*W_op^2), 1, 1, 1, 1/5, 1/5, 1/5];
params.gamma2 = ones(1,9)/8;

disp('Algorithm starts.')
if params.HR_sparse == 0 && params.HR_stripe_rate == 0
    params.gamma1 = [1/(2+32*W_op^2) 1/(1+32*W_op^2)];
    params.gamma2 = ones(1,6)/2;
    STfusion;
elseif params.HR_sparse ~= 0 && params.HR_stripe_rate == 0
    params.gamma1 = [1/(2+32*W_op^2) 1/(1+32*W_op^2) 1 1 1];
    params.gamma2 = ones(1,6)/5;
    STfusion_with_sparsenoise;
elseif params.HR_sparse == 0 && params.HR_stripe_rate ~= 0
    params.gamma1 = [1/(2+32*W_op^2) 1/(1+32*W_op^2) 1/5 1/5 1/5];
    params.gamma2 = ones(1,9)/5;
    STfusion_with_stripenoise;
elseif params.HR_sparse ~= 0 && params.HR_stripe_rate ~= 0
    params.gamma1 = [1/(2+32*W_op^2), 1/(1+32*W_op^2), 1, 1, 1, 1/5, 1/5, 1/5];
    params.gamma2 = ones(1,9)/8;
    STfusion_with_sparsenoise_and_stripenoise;
end