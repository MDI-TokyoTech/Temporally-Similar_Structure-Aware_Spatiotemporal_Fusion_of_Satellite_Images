% Xは4DテンソルでX(y,x,b,d)には、y座標、x座標、bバンド、方向dの差分値が入っているとする
% dirはl2ノルムをまとめる方向
% dir = 3ならバンド方向にだけまとめる
% dir = [3 4]ならバンド方向と差分方向合わせてまとめる

function result = prox12band(X, gamma, dir)
    if nargin < 3
        dir = [3,4];
    end
    l2X = sqrt(sum(X.^2, dir));
    T = max(1 - gamma./l2X, 0);
    result = T.*X;
end