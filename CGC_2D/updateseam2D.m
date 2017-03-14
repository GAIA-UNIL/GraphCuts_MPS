function [seamnode,costmap] = updateseam2D(seamnode,newseamnode,costmap,Nadd,lineno,check)

% This function is to update the seam nodes. The process consist:
% 1. remove the old list; 2. check the cut label; 3. add the new seam node
% to seamnode; 4 update costmap
% seamnode have 10 columns.
% x,y,z,cutcost,ind0,ind1,TIind0,TIind0adjacent,TIind1,TIind1adjacent

if nargin < 4
    Nadd = 0;
    lineno = [];
    check = nan(0,7);
end

if Nadd > 0
    % if no previous cut exists in the overlap area. just add the new seam
    % nodes
    checkseam = nan(0,10);
    seamnode_old = seamnode(lineno,:);
    seamnode(lineno,:) = [];
    costmap(seamnode_old(:,5:6)) = 0;

    for nn = 1:Nadd
        %% check the arc of the seam node
        ch = [num2str(check(nn,1)),num2str(check(nn,2)),num2str(check(nn,3))];
%        disp(ch);
        switch ch
            case '000'
                % keep the old cut which means find the corronsponding line
                % from the seamnode_old and add this line to the new seamnode
                % since it is remove in line 8
                checkseam = [checkseam;seamnode_old(nn,:)];
            case '110'
                % make a new cut at the same location with the value of the
                % cost of patch B and C, patch C take place patch A, retain
                % patch B
                changeline = find(newseamnode(:,5)==seamnode_old(nn,6) & newseamnode(:,6)==seamnode_old(nn,5));
                newseamnode(changeline,4) = check(nn,7);
                % exturn the position of patch B since the cut label extruned
                newseamnode(changeline,7) = seamnode_old(nn,10);
                newseamnode(changeline,8) = seamnode_old(nn,9);
                % The value of the old patch keep the same
            case '101'
                % make a new cut at the same location with the value of the
                % cost of patch A and patch C, patch C take place of patch
                % B, retain patch A
                changeline = find(newseamnode(:,5)==seamnode_old(nn,5) & newseamnode(:,6)==seamnode_old(nn,6));
                newseamnode(changeline,4) = check(nn,6);
            case '111'
                % the old seam is removed by the new patch
            otherwise
                disp('unexpected cut, check!');
                pause;
        end % different cases
    end % Loop for all the added seam nodes
    
    %% update the seamnode
    seamnode = [seamnode;newseamnode;checkseam];
    costmap(checkseam(:,5)) = 0.5 * checkseam(:,4);
    costmap(checkseam(:,6)) = 0.5 * checkseam(:,4);
%    abs(checkseam(:,7)-checkseam(:,9));
%    costmap(checkseam(:,6)) = abs(checkseam(:,8)-checkseam(:,10)) + eps;
else
    seamnode = [seamnode;newseamnode];    
end
costmap(newseamnode(:,5)) = 0.5 * newseamnode(:,4);
costmap(newseamnode(:,6)) = 0.5 * newseamnode(:,4);
% costmap(newseamnode(:,5)) = abs(newseamnode(:,7)-newseamnode(:,9)) + eps;
% costmap(newseamnode(:,6)) = abs(newseamnode(:,8)-newseamnode(:,10)) + eps;
end