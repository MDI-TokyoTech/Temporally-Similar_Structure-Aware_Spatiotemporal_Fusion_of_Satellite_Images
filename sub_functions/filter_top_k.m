function A_filtered = filter_top_k(A, k)
    % A: 任意サイズの (H × W × 4) の3D配列
    % k: 残す要素数（1, 2, 3, 4 のいずれか）

    % 配列サイズを取得
    [H, W, D] = size(A);

    % k が 1, 2, 3, 4 以外ならエラー
    if ~ismember(k, [1, 2, 3, 4])
        error('k must be 1, 2, 3, or 4');
    end

    % k = 4 の場合、A をそのまま返す
    if k == D
        A_filtered = A;
        return;
    end

    % 各 (y, x) の 4 要素を降順ソートし、最大 k 個のインデックスを取得
    [~, sortedIdx] = sort(A, 3, 'descend');

    % 上位 k 個のインデックス
    idx_top_k = sortedIdx(:,:,1:k);

    % マスクを作成 (すべて false で初期化)
    mask = false(H, W, D);

    % 各次元のインデックスを作成
    rowIdx = repmat((1:H)', [1, W, k]); % y インデックス (H×W×k)
    colIdx = repmat(reshape(1:W, 1, W), [H, 1, k]); % x インデックス (H×W×k)

    % 3次元の線形インデックスを取得し、mask を true に設定
    mask(sub2ind(size(A), rowIdx, colIdx, idx_top_k)) = true;

    % マスクを適用して、それ以外の要素を 0 にする
    A_filtered = A .* mask;
end
