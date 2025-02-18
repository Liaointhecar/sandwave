function index = new_find(R)
index = [];
index(1) = find(R==max(max(R)));
search_end = 0;
search_start = 0;
s_r = size(R);

% i 撞到邊界則停止
boundary_up = 1:size(R,1):numel(R)-size(R,1)+1;
boundary_down = size(R,1):size(R,1):numel(R);
boundary_left = 1:size(R,1);
boundary_right = numel(R)-size(R,1)+1:numel(R);

% go down 
while search_start <= search_end
    search_start = search_end+1;
    search_end = length(index); % 一直找到最新一個
    for g = search_start:search_end % index 編號
        i = index(g); % index 矩陣上的位置
        % down
        if ismember(i,boundary_down) == 0 % 本身非邊界
            if R(i+1) > 0 && ismember(i+1,index) == 0 % R(i) 矩陣上的值
                % 下一個大於零、不是已有的index之一
                s = length(index)+1; % index 目前有幾個
                index(s) = i+1; % 加上新的index
            end
        end
        % left
        if ismember(i,boundary_left) == 0 % 本身非邊界
            if  R(i-s_r(2)) > 0 && ismember(i-s_r(2),index) == 0 % R(i) 矩陣上的值
                % 下一個大於零、不是已有的index之一
                s = length(index)+1; % index 目前有幾個
                index(s) = i-s_r(2); % 加上新的index
            end
        end
        % up
        if ismember(i,boundary_up) == 0 % 本身非邊界
            if R(i-1) > 0 && ismember(i-1,index) == 0  % R(i) 矩陣上的值
                % 下一個大於零、不是已有的index之一
                s = length(index)+1; % index 目前有幾個
                index(s) = i-1; % 加上新的index
            end
        end
        % right
        if ismember(i,boundary_right) == 0 % 本身非邊界
            if R(i+s_r(2)) > 0 && ismember(i+s_r(2),index) == 0  % R(i) 矩陣上的值
                % 下一個大於零、不是已有的index之一
                s = length(index)+1; % index 目前有幾個
                index(s) = i+s_r(2); % 加上新的index
            end
        end
    end
end
