function index = new_find(R)
index = [];
index(1) = find(R==max(max(R)));
search_end = 0;
search_start = 0;
s_r = size(R);

% i ������ɫh����
boundary_up = 1:size(R,1):numel(R)-size(R,1)+1;
boundary_down = size(R,1):size(R,1):numel(R);
boundary_left = 1:size(R,1);
boundary_right = numel(R)-size(R,1)+1:numel(R);

% go down 
while search_start <= search_end
    search_start = search_end+1;
    search_end = length(index); % �@�����̷s�@��
    for g = search_start:search_end % index �s��
        i = index(g); % index �x�}�W����m
        % down
        if ismember(i,boundary_down) == 0 % �����D���
            if R(i+1) > 0 && ismember(i+1,index) == 0 % R(i) �x�}�W����
                % �U�@�Ӥj��s�B���O�w����index���@
                s = length(index)+1; % index �ثe���X��
                index(s) = i+1; % �[�W�s��index
            end
        end
        % left
        if ismember(i,boundary_left) == 0 % �����D���
            if  R(i-s_r(2)) > 0 && ismember(i-s_r(2),index) == 0 % R(i) �x�}�W����
                % �U�@�Ӥj��s�B���O�w����index���@
                s = length(index)+1; % index �ثe���X��
                index(s) = i-s_r(2); % �[�W�s��index
            end
        end
        % up
        if ismember(i,boundary_up) == 0 % �����D���
            if R(i-1) > 0 && ismember(i-1,index) == 0  % R(i) �x�}�W����
                % �U�@�Ӥj��s�B���O�w����index���@
                s = length(index)+1; % index �ثe���X��
                index(s) = i-1; % �[�W�s��index
            end
        end
        % right
        if ismember(i,boundary_right) == 0 % �����D���
            if R(i+s_r(2)) > 0 && ismember(i+s_r(2),index) == 0  % R(i) �x�}�W����
                % �U�@�Ӥj��s�B���O�w����index���@
                s = length(index)+1; % index �ثe���X��
                index(s) = i+s_r(2); % �[�W�s��index
            end
        end
    end
end
