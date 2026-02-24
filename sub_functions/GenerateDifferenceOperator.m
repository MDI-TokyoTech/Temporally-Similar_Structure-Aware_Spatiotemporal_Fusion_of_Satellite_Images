function [D, Dt] = GenerateDifferenceOperator()

    up = @(z) [z(2:end,:,:); z(1,:,:)];
    down = @(z) [z(end,:,:); z(1:end-1,:,:)];
    right = @(z) [z(:,end,:) z(:,1:end-1,:)];
    left = @(z) [z(:,2:end,:) z(:,1,:)];
    
%     cellD{1} = @(z) down(right(z)) -z;
%     cellD{2} = @(z) down(z) -z;
%     cellD{3} = @(z) down(left(z)) -z;
%     cellD{4} = @(z) right(z) -z;
%     cellD{5} = @(z) left(z) -z;
%     cellD{6} = @(z) up(right(z)) -z;
%     cellD{7} = @(z) up(z) -z;
%     cellD{8} = @(z) up(left(z)) -z;
%     
%     cellDt{1} = @(z) cellD{8}(z);
%     cellDt{2} = @(z) cellD{7}(z);
%     cellDt{3} = @(z) cellD{6}(z);
%     cellDt{4} = @(z) cellD{5}(z);
%     cellDt{5} = @(z) cellD{4}(z);
%     cellDt{6} = @(z) cellD{3}(z);
%     cellDt{7} = @(z) cellD{2}(z);
%     cellDt{8} = @(z) cellD{1}(z);

    Dv = @(z) up(z) -z;
    Dh = @(z) left(z) -z;
    Dbr = @(z) up(left(z)) -z;
    Dbl = @(z) up(right(z)) -z;

    Dvt = @(z) down(z) -z;
    Dht = @(z) right(z) -z;
    Dbrt = @(z) down(right(z)) -z;
    Dblt = @(z) down(left(z)) -z;

    D = @(z) cat(4, Dv(z), Dh(z), Dbr(z), Dbl(z));
    Dt = @(z) Dvt(z(:,:,:,1)) + Dht(z(:,:,:,2)) + Dbrt(z(:,:,:,3)) + Dblt(z(:,:,:,4));

end

% 1 2 3
% 4 * 5
% 6 7 8

% 
% function[y] = slide_upleft(x)
% %     up_x = x(1,2:end,:);
% %     left_x = circshift(x(:,1,:), [-1, 0]);
% %     x = [[x(2:end, 2:end,:); up_x] left_x];
%     y = [[x(2:end, 2:end,:); x(1,2:end,:)] circshift(x(:,1,:), [-1, 0])];
% end
% 
% function[y] = slide_upright(x)
%     y = [circshift(x(:,end,:), [-1, 0]) [x(2:end, 1:end-1,:); x(1,1:end-1,:)]];
% end
% 
% function[y] = slide_downleft(x)
%     y = [circshift(x(:,end,:), [1, 0]) [x(end, 1:end-1,:); x(1:end-1,1:end-1,:)]];
% end
% 
% function[y] = slide_downright(x)
%     y = [[x(end, 2:end,:); x(1:end-1, 2:end,:)] circshift(x(:,1,:), [1, 0])];
% end