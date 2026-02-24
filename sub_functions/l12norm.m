function result = l12norm(X, dir)
    if nargin < 2
        l2X = sqrt(sum(X.^2, [3, 4]));
    else
        l2X = sqrt(sum(X.^2,dir));
    end
    result = sum(l2X(:));
end

