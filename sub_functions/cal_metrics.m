function [p,r,ss,cc,sm] = cal_metrics(pre, gt)

    % --- 画像サイズ確認 ---
    [H, W, C] = size(gt);

    % --- 値が有効なマスク（各バンド単位）---
    valid_mask_band = (gt >= 0) & (gt <= 1);  % size: H x W x C

    % --- SAM用：ピクセル単位で全バンド有効なマスク ---
    valid_mask_pixel = all(valid_mask_band, 3);  % size: H x W

    % =======================
    % === PSNR / RMSE / CC ===
    % =======================
    p_list = zeros(1, C);
    rmse_list = zeros(1, C);
    pre_all = []; gt_all = [];

    for c = 1:C
        mask = valid_mask_band(:,:,c);
        pre_c = pre(:,:,c);
        gt_c  = gt(:,:,c);

        pre_valid = pre_c(mask);
        gt_valid  = gt_c(mask);

        if ~isempty(gt_valid)
            p_list(c) = psnr(pre_valid, gt_valid);
            rmse_list(c) = sqrt(mean((pre_valid - gt_valid).^2));
            pre_all = [pre_all; pre_valid(:)];
            gt_all  = [gt_all;  gt_valid(:)];
        else
            p_list(c) = NaN;
            rmse_list(c) = NaN;
        end
    end

    p = mean(p_list, 'omitnan');
    r = mean(rmse_list, 'omitnan');

    if isempty(gt_all)
        cc = NaN;
    else
        C_corr = corrcoef(pre_all, gt_all);
        cc = C_corr(1,2);
    end

    % ============
    % === SSIM ===
    % ============
    ss_list = zeros(1, C);
    for c = 1:C
        pre_c = pre(:,:,c);
        gt_c  = gt(:,:,c);

        % 無効なインデックス（範囲外）
        invalid_mask = (gt_c < 0) | (gt_c > 1);

        % gtの範囲外値をpreの値で置換
        gt_filled = gt_c;
        gt_filled(invalid_mask) = pre_c(invalid_mask);

        % SSIM計算
        ss_list(c) = ssim(pre_c, gt_filled);
    end
    ss = mean(ss_list, 'omitnan');

    % ===============
    % === SAM計算 ===
    % ===============
    % 有効ピクセルだけ抽出
    pre_valid = pre(repmat(valid_mask_pixel, [1 1 C]));
    gt_valid  = gt(repmat(valid_mask_pixel, [1 1 C]));

    pre_valid = reshape(pre_valid, [], C);
    gt_valid  = reshape(gt_valid, [], C);

    dot_prod = sum(pre_valid .* gt_valid, 2);
    norm_pre = sqrt(sum(pre_valid.^2, 2));
    norm_gt  = sqrt(sum(gt_valid.^2, 2));
    sm = mean(acos(dot_prod ./ (norm_pre .* norm_gt)), 'omitnan');
end

% === 無効値を4近傍の平均で補間する関数 ===
function filled = fill_invalid_with_neighbor_avg(img, valid_mask)
    filled = img;
    [H, W] = size(img);

    invalid_idx = find(~valid_mask);

    for idx = invalid_idx'
        [i, j] = ind2sub([H, W], idx);

        neighbors = [];

        if i > 1     && valid_mask(i-1,j), neighbors(end+1) = img(i-1,j); end
        if i < H     && valid_mask(i+1,j), neighbors(end+1) = img(i+1,j); end
        if j > 1     && valid_mask(i,j-1), neighbors(end+1) = img(i,j-1); end
        if j < W     && valid_mask(i,j+1), neighbors(end+1) = img(i,j+1); end

        if ~isempty(neighbors)
            filled(i,j) = mean(neighbors);
        else
            filled(i,j) = 0; % 近傍全て無効な場合は0で埋める（任意対応）
        end
    end
end
