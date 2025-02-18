%%Based on the method developed by Duffy and Hughes Clarke (2005)
size_of_vq2016 = size(vq_2016);
size_of_vq2017 = size(vq_2017);
% get_min = min(min(size_of_vq2016),min(size_of_vq2017));
get_max = max(max(size_of_vq2016),max(size_of_vq2017));

magical_number = 1;



%%
window_size = 220; % set window size in meter (small)
w = window_size/sr;


%set start
%x1 = get_max - magical_number - w + 1;
x1 = 1950;
y1 = 50; %fix

%set small(t1) and big(t2) window
range = x1 - y1; 

% for cross-correlation
A1_2016 = vq_2016(x1-range:x1,y1:y1+range); % small size
A1_2017 = vq_2017(x1-range-(w-1):x1+(w-1),...
    y1-(w-1):y1+range+(w-1)); % big size

%%%%for plot%%%%
B1_2016 = f_2016(x1-range:x1,y1:y1+range);
B1_2017 = f_2017(x1-range:x1,y1:y1+range);

a = A1_2016;
b = A1_2017;

w_s = 1:w; % 1st (small window)
w_b = 1:3*w-2; % 2nd (big window)
s_all = size(a);

%% 
clear r R r_full2...
    displacement_u1 displacement_u2 displacement_u3...
    displacement_v1 displacement_v2 displacement_v3...
    distance_1 distance_2 distance_3 distance_4
displacement_u1(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
displacement_v1(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
displacement_u2(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
displacement_v2(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
displacement_u3(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
displacement_v3(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
R_matrix_440(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;


azimuth_0(1:floor(s_all(1)/(w/2))-1,1:floor(s_all(2)/(w/2))-1) = 0;
for i = 1:floor(s_all(1)/(w/2))-1
    for j = 1:floor(s_all(2)/(w/2))-1
  
        f_corr2 = a(w_s + floor((i-1)*(w/2)), w_s + floor((j-1)*(w/2)));
        g_corr2 = b(w_b + floor((i-1)*(w/2)), w_b + floor((j-1)*(w/2))); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp_f_corr = f_corr2(:)';
%         if var(tmp_f_corr) < 0.35 %0.35
        if var(tmp_f_corr) < 0.35 || max(tmp_f_corr==0)>0 %0.35    
            displacement_u1(i,j) = 0;
            displacement_v2(i,j) = 0;
            displacement_u2(i,j) = 0;
            displacement_v2(i,j) = 0;
            displacement_u3(i,j) = 0;
            displacement_v3(i,j) = 0;
        else
            r = normxcorr2(f_corr2,g_corr2);
            R = r(w-1:length(r)-w+2,w-1:length(r)-w+2); % cut edge
        
            %R = r(w-1-5:length(r)-w+2+5,w-1-5:length(r)-w+2+5) % cut edge
            %R = r  %take all
            
            %Threshold~~~
            %R(R<0.6) = 0
            R(R<1/sqrt(2)) = 0;
            %R(R<0.9) = 0

  
            % Method i --- Maximum correlation %
             
            [r1, c1] = find(R==max(max(R)));
            if length(r1)>1 || length(c1)>1
                displacement_u1(i,j) = 0;
                displacement_v1(i,j) = 0;
                displacement_u2(i,j) = 0;
                displacement_v2(i,j) = 0;
                displacement_u3(i,j) = 0;
                displacement_v3(i,j) = 0;
                R_matrix_440(i,j) = 0;
            else
                displacement_u1(i,j) = c1-0.5*(length(R)+1);
                displacement_v1(i,j) = r1-0.5*(length(R)+1);
                R_matrix_440(i,j) = max(max(R));
        
                clear r1 c1
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                index = new_find(R);
                [h_r, h_c] = ind2sub(size(R),index);
                h_r = h_r';
                h_c = h_c';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % [h_r h_c] = find(R>0)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Method ii --- Weighted centroid %
                %h_r = index(:,1);
                %h_c = index(:,2);
                cen = (length(R)+1)/2;
                h_r1 = h_r-cen;
                h_c1 = h_c-cen;
                i_h = 1:length(h_c);
                h_x_s = sum(h_c1(i_h).*diag(R(h_r,h_c)))./ sum(diag(R(h_r,h_c)));
                h_y_s = sum(h_r1(i_h).*diag(R(h_r,h_c)))./ sum(diag(R(h_r,h_c)));
                h_x_s(isnan(h_x_s)) = 0;
                h_y_s(isnan(h_y_s)) = 0;
                
                displacement_u2(i,j) = h_x_s;
                displacement_v2(i,j) = h_y_s;
                
                % Method iii --- Regression line %
                X = [ones(length(h_c1),1)  h_c1];
                reg = X\h_r1;
                %----------------------------------------------------------------start
                displacement_u3(i,j) = (-1).*reg(1)./(reg(2)+(1./reg(2)));

%                 if reg(2) == 0
%                     displacement_v3(i,j) = 0; 
%                 else
                    displacement_v3(i,j) = (-1).*displacement_u3(i,j)./reg(2);
%                     d_r_y2 = reg(1)+ reg(2).*((max(h_c1)+min(h_c1))./2);
%                     displacement_u3(i,j) = (max(h_c1)+min(h_c1))./2;
%                     displacement_v3(i,j) = d_r_y2;
%                     azimuth_0(i,j) = reg(2);
%                 end
                %----------------------------------------------------------------end
                %r_full2(1+(2*w-1).*(i-1):1+(2*w-1).*(i-1)+length(R)-1,1+(2*w-1).*(j-1):1+(2*w-1).*(j-1)+length(R)-1) = R
            end
        end
    end
end

%migration distance
distance_1 = sqrt(displacement_u1.^2+displacement_v1.^2);
distance_2 = sqrt(displacement_u2.^2+displacement_v2.^2);
distance_3 = sqrt(displacement_u3.^2+displacement_v3.^2);



%bias
displacement_u1(distance_1>w) = 0;
displacement_v1(distance_1>w) = 0;

displacement_u2(distance_2>w) = 0;
displacement_v2(distance_2>w) = 0;

displacement_u3(distance_3>w) = 0;
displacement_v3(distance_3>w) = 0;


std_440 = [std(displacement_u2(displacement_u2>0))*5 std(displacement_v2(displacement_v2>0))*5 std(distance_2(distance_2>0))*5];
%%
std_x=[55 std_55(1);80 std_80(1);110 std_110(1);150 std_150(1);170 std_170(1);180 std_180(1);200 std_200(1);220 std_220(1);230 std_230(1);250 std_250(1);330 std_330(1);440 std_440(1)];
std_y=[55 std_55(2);80 std_80(2);110 std_110(2);150 std_150(2);170 std_170(2);180 std_180(2);200 std_200(2);220 std_220(2);230 std_230(2);250 std_250(2);330 std_330(2);440 std_440(2)];
std_d=[55 std_55(3);80 std_80(3);110 std_110(3);150 std_150(3);170 std_170(3);180 std_180(3);200 std_200(3);220 std_220(3);230 std_230(3);250 std_250(3);330 std_330(3);440 std_440(3)];
std_9=[55 (1-sum(R_matrix_55(:)>=0.9)/sum(R_matrix_55(:)>0))*100;80 (1-sum(R_matrix_80(:)>=0.9)/sum(R_matrix_80(:)>0))*100;110 (1-sum(R_matrix_110(:)>=0.9)/sum(R_matrix_110(:)>0))*100;150 (1-sum(R_matrix_150(:)>=0.9)/sum(R_matrix_150(:)>0))*100;170 (1-sum(R_matrix_170(:)>=0.9)/sum(R_matrix_170(:)>0))*100;180 (1-sum(R_matrix_180(:)>=0.9)/sum(R_matrix_180(:)>0))*100;200 (1-sum(R_matrix_200(:)>=0.9)/sum(R_matrix_200(:)>0))*100;220 (1-sum(R_matrix_220(:)>=0.9)/sum(R_matrix_220(:)>0))*100;230 (1-sum(R_matrix_230(:)>=0.9)/sum(R_matrix_230(:)>0))*100;250 (1-sum(R_matrix_250(:)>=0.9)/sum(R_matrix_250(:)>0))*100;330 (1-sum(R_matrix_330(:)>=0.9)/sum(R_matrix_330(:)>0))*100;440 (1-sum(R_matrix_440(:)>=0.9)/sum(R_matrix_440(:)>0))*100];

figure;hold on;yyaxis left;plot(std_x(:,1),std_x(:,2));plot(std_y(:,1),std_y(:,2));yyaxis right;plot(std_9(:,1),std_9(:,2));

%% plot %
method = '3';
displacement_u = strcat('displacement_u',method);
displacement_v = strcat('displacement_v',method);

[xq_A1,yq_A1] = meshgrid(1:length(a),1:length(a));

i = figure;
axis tight manual % this ensures that getframe() returns a consistent size
mesh(xq_A1,yq_A1,B1_2016);
hold on;
x_component = floor(w/2):floor(w/2):floor(w/2)*(floor(s_all(1)./floor(length(w_s)./2))-1);
[vector_x,vector_y] = meshgrid(x_component,x_component);
hold on;
quiver(vector_x,vector_y,eval(displacement_u),eval(displacement_v),'color','black')
set(gca,'FontSize',9)
hold off
axis equal
az = 0;
el = 90;
view(az, el);
title('m2 2016 black','fontsize',9)
xticks(0:250:2500);
xticklabels((xticks+50)*5+min(b_2016(:,1)));
yticks(0:250:2500);
yticklabels((yticks+50)*5+min(b_2016(:,2)));
%set(h,'Units','centimeters','position',[5,5,10,20])
   drawnow 
      % Capture the plot as an image 
      frame = getframe(i); 
      im = frame2im(frame); 
      [~,~] = rgb2ind(im,256);

%       image1_name = 'm2_2016_black.png'
%       image1_path = strcat(area_path,'\',image1_name);
%       imwrite(imind,cm,image1_path);
%~~~~~finish~~~~~%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method2 separate by orientation %%

count2 = size(displacement_u2);
displacement_u2_1 = zeros(count2(1),count2(2));
displacement_v2_1 = zeros(count2(1),count2(2));
displacement_u2_2 = zeros(count2(1),count2(2));
displacement_v2_2 = zeros(count2(1),count2(2));
displacement_u2_3 = zeros(count2(1),count2(2));
displacement_v2_3 = zeros(count2(1),count2(2));
displacement_u2_4 = zeros(count2(1),count2(2));
displacement_v2_4 = zeros(count2(1),count2(2));

direction = zeros(count2(1),count2(2));

for i = 1:count2(1)
    for j = 1:count2(2)
    if displacement_u2(i,j)*displacement_v2(i,j) < 0
        if displacement_u2(i,j) < 0
            displacement_u2_2(i,j) = displacement_u2(i,j); %red
            displacement_v2_2(i,j) = displacement_v2(i,j);
            direction(i,j) = (atan(displacement_u2(i,j)/displacement_v2(i,j)))*180./pi; % direction
        else 
            displacement_u2_4(i,j) = displacement_u2(i,j); %white
            displacement_v2_4(i,j) = displacement_v2(i,j);           
            direction(i,j) = (atan(displacement_u2(i,j)/displacement_v2(i,j)))*180./pi;% + 180;
        end
    else
        if displacement_u2(i,j) > 0
            displacement_u2_1(i,j) = displacement_u2(i,j); %black
            displacement_v2_1(i,j) = displacement_v2(i,j);
            direction(i,j) = (atan(displacement_u2(i,j)/displacement_v2(i,j)))*180./pi;
        else
            displacement_u2_3(i,j) = displacement_u2(i,j); %yellow
            displacement_v2_3(i,j) = displacement_v2(i,j);
            direction(i,j) = (atan(displacement_u2(i,j)/displacement_v2(i,j)))*180./pi ;% 180;
        end
    end
    end
end

[xq_A1,yq_A1] = meshgrid(1:length(a),1:length(a));
%%
count2 = size(displacement_u2);
direction3A = zeros(count2(1),count2(2));
direction3A_abs = zeros(count2(1),count2(2));
for i = 1:count2(1)
    for j = 1:count2(2)
direction3A_abs(i,j) = sqrt(displacement_u3(i,j)^2 + displacement_v3(i,j)^2);
direction3A(i,j) = 90-atan2(displacement_u3(i,j)/direction3A_abs(i,j),displacement_v3(i,j)/direction3A_abs(i,j))*180/pi;
    end
end
dir_dis=[direction3A(:) distance_3(:)];
%% Method3 separate by orientation %%

count2 = size(displacement_u3);
displacement_u3_1 = zeros(count2(1),count2(2));
displacement_v3_1 = zeros(count2(1),count2(2));
displacement_u3_2 = zeros(count2(1),count2(2));
displacement_v3_2 = zeros(count2(1),count2(2));
displacement_u3_3 = zeros(count2(1),count2(2));
displacement_v3_3 = zeros(count2(1),count2(2));
displacement_u3_4 = zeros(count2(1),count2(2));
displacement_v3_4 = zeros(count2(1),count2(2));

direction3 = zeros(count2(1),count2(2));

for i = 1:count2(1)
    for j = 1:count2(2)
    if displacement_u3(i,j)*displacement_v3(i,j) < 0
        if displacement_u3(i,j) < 0
            displacement_u3_2(i,j) = displacement_u3(i,j); %red
            displacement_v3_2(i,j) = displacement_v3(i,j);
            direction3(i,j) = (atan(displacement_u3(i,j)/displacement_v3(i,j)))*180./pi; % direction
        else 
            displacement_u3_4(i,j) = displacement_u3(i,j); %white
            displacement_v3_4(i,j) = displacement_v3(i,j);           
            direction3(i,j) = (atan(displacement_u3(i,j)/displacement_v3(i,j)))*180./pi+ 180;
        end
    else
        if displacement_u3(i,j) > 0
            displacement_u3_1(i,j) = displacement_u3(i,j); %black
            displacement_v3_1(i,j) = displacement_v3(i,j);
            direction3(i,j) = (atan(displacement_u3(i,j)/displacement_v3(i,j)))*180./pi;
        else
            displacement_u3_3(i,j) = displacement_u3(i,j); %yellow
            displacement_v3_3(i,j) = displacement_v3(i,j);
            direction3(i,j) = (atan(displacement_u3(i,j)/displacement_v3(i,j)))*180./pi+ 180;
        end
    end
    end
end

[xq_A1,yq_A1] = meshgrid(1:length(a),1:length(a));
%% plot separate by orientation%
i = figure;
axis tight manual % this ensures that getframe() returns a consistent size
mesh(xq_A1,yq_A1,B1_2016);
hold on



[vector_x,vector_y] = meshgrid(x_component,x_component);
L = quiver(vector_x,vector_y,displacement_u3_1,displacement_v3_1,'color','black');
quiver(vector_x,vector_y,displacement_u3_2,displacement_v3_2,'color','red');
quiver(vector_x,vector_y,displacement_u3_3,displacement_v3_3,'color','yellow');
quiver(vector_x,vector_y,displacement_u3_4,displacement_v3_4,'color','white');
set(gca,'FontSize',18)
hold off
axis equal
az = 0;
el = 90;
view(az, el);
title('m2 2016','fontsize',18)
%set(h,'Units','centimeters','position',[5,5,10,20])
   drawnow 
      % Capture the plot as an image 
      frame = getframe(i); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256);

%       image2_name = 'm2_2016.png'
%       image2_path = strcat(area_path,'\',image2_name);
%       imwrite(imind,cm,image2_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
distance_3_c = distance_3;
distance_3_c(isnan(direction3A))=nan;
distance_3_c = distance_3_c(~isnan(distance_3_c));
direction_3_c = direction3A(~isnan(direction3A));
dir_dis_rose = [direction_3_c distance_3_c.*sr];
dir_dis_rose_100 = dir_dis_rose(~any(dir_dis_rose(:,2)>100, 2), :);
dir_dis_rose_5_100 = dir_dis_rose_100(~any(dir_dis_rose_100(:,2)<1, 2), :);


[figure_handle, count, speeds, directions, Table] = WindRose(dir_dis_rose_5_100(:,1)-30, dir_dis_rose_5_100(:,2), 'cMap', 'jet', 'FreqLabelAngle', 'auto');
dis_5_N = dir_dis_rose_5_100(dir_dis_rose_5_100(:,1)<90|dir_dis_rose_5_100(:,1)>=270,:);
dis_5_S = dir_dis_rose_5_100(dir_dis_rose_5_100(:,1)>=90&dir_dis_rose_5_100(:,1)<270,:);
dis_5_all = [dis_5_N(:,2);dis_5_S(:,2)*(-1)];
figure;histogram(dis_5_all,40);

%%%%  b i g  %%%%
%% back to TM2 %%
clear arrow_distance
%arrow_distance = [(vector_x(:)+50)*5+min(b_2016(:,1)), (vector_y(:)+50)*5+min(b_2016(:,2)), distance_2(:).*sr];
% arrow_distance = [(vector_x(:)+50)*5+super_min_x, (vector_y(:)+50)*5+super_min_y, direction(:), distance_2(:).*sr];
% arrow_distance = [vector_x(:)*sr+super_min_x, vector_y(:)*sr+super_min_y, vector_x(:), vector_y(:), direction(:), distance_2(:).*sr];



%filename = strcat(area_name,'_distance.xlsx');
%filename_path = strcat(area_path,'\',filename);
% writematrix(arrow_distance,filename_path,'Sheet',1,'Range','A1');
