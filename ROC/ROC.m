fid = fopen('breastinfo_simple.txt','r');
breastID = str2double(fgetl(fid));
s1 = str2double(fgetl(fid));
s2 = str2double(fgetl(fid));
s3 = str2double(fgetl(fid));
class = str2double(fgetl(fid));
fclose(fid);

load mtype.mat;
load pval.mat;

muscle_wall = 153;
skin_start = 138;   

% Convert vector into cube
mtype_cube = zeros(s1,s2,s3); % each voxel is .5mmx.5mmx.5mm
pval_cube = zeros(s1,s2,s3);
cur_pos = 1;
for k=1:s3
    for j=1:s2
        for i= 1:s1
            mtype_cube(i,j,k) = mtype(cur_pos);
            pval_cube(i,j,k) = pval(cur_pos);
            cur_pos = cur_pos + 1;
        end 
    end
end

% subsample cubes in order to solve sparse matrix
s1_ss = floor(s1/2); % voxels are now 1mmx1mmx1mm
s2_ss = floor(s2/2);
s3_ss = floor(s3/2);
xi = 1; yi = 1; zi = 1;
mtype_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
pval_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
for z=1:2:s3-1
    for y = 1:2:s2
        for x = 1:2:s1
            mid = mtype_cube(x,y,z);
            pid = pval_cube(x,y,z);
            mtype_cube_subsamp(xi,yi,zi) = mid;
            pval_cube_subsamp(xi,yi,zi) = pid;
            xi = xi+1;
        end
        xi = 1;
        yi = yi + 1;
    end
    yi = 1;
    zi = zi + 1;
end
% some voxels not converted to muscle during subsampling
% so do that now
% still need to figure out how to get pval converted for fdtd
for x=1:s1_ss
    for y=1:s2_ss
        for z=1:s3_ss
            if x > 153
               mtype_cube_subsamp(x,y,z) = -4;
                pval_cube_subsamp(x,y,z) = 1;
            end
        end
    end
end

% save mtype_cube_subsamp; save pval_cube_subsamp.mat;
% load mtype_cube_subsamp.mat;
model = mtype_cube_subsamp; air_id = -1;
[s1_ss, s2_ss, s3_ss] = size(model);
%model = model(18:end,:,:);
%[s1_ss, s2_ss, s3_ss] = size(model);


% file_name = 'scat_fib_1.csv';
% create_csv(file_name,model,s1_ss,s2_ss,s3_ss,air_id,mtype_cube_subsamp,pval_cube_subsamp);

figure; colormap(gray);
contourf(model(:,:,floor(s3_ss/2))); 

% model_pval = pval_cube_subsamp;
figure;
colormap jet;
contourf((0:s2_ss-1)*.1,(0:s1_ss-1)*.1,model(:,:,floor(s3_ss/2)),'LineStyle','none');
% colormap(flipud(gray));brighten(.4);
xlabel('Distance (cm)','FontSize',14); ylabel('Distance (cm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);

% Create Norminal temperature model
tumor_on = 0; tumor_depth = 8;  tumor_radius = 10; Tambient = 27; Tart = 37; 
% muscle_wall = 154; skin_start = 0; % can't remember why I put this in
[T_3d_nom,tissue_3d_nom] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start);
% Generate A_norminal
tum_x_cen =0; tum_y_cen =0; tum_z_cen = 0;
A_norminal = gen_A(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);
A_index = zeros(1,17*176*153);
i = 1;
for z = 1:s3_ss
    for x = 1:17
        for y = 1:s2_ss
            A_index(1,i) = (z-1)*158*176+(x-1)*176+y;
            i = i+1;
        end
    end
end
A_norminal(:,A_index)=[];
A_norminal(A_index',:)=[];

figure;
colormap jet;
contourf(T_3d_nom(:,:,floor(s3_ss/2)));
colorbar;
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
title('Normal Temperature Profile (\circC)','FontSize',14);

% Generate Anomaly temperatures with tumor radius = 10 mm, tumor depth = 10 mm
tumor_on = 1; tum_y_cen = 90; tum_z_cen = floor(s3_ss/2);

tumor_radius = 10; tumor_depth = 10; tum_x_cen = 20 + tumor_depth;% tumor dept is 1cm
[T_3d_abn1,tissue_3d_abn1] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);

tumor_radius = 10; tumor_depth = 20; tum_x_cen = 20 + tumor_depth;% tumor dept is 2cm
[T_3d_abn2,tissue_3d_abn2] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);

figure; 
colormap jet;
contourf(T_3d_abn1(:,:,floor(s3_ss/2)));
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);
title('Temperature Profile of 1cm Deep Tumor (\circC)','FontSize',14);

figure; 
colormap jet;
contourf(T_3d_abn2(:,:,floor(s3_ss/2)));
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);
title('Temperature Profile of 2cm Deep Tumor (\circC)','FontSize',14);

for z = 1:s3_ss
    for x = 18:s1_ss
        for y = 1:s2_ss
            T_3d_nom(x-17,y,z) = T_3d_nom(x,y,z);
            T_3d_abn1(x-17,y,z) = T_3d_abn1(x,y,z);
            T_3d_abn2(x-17,y,z) = T_3d_abn2(x,y,z);
            
        end
    end
end

model = model(18:end,:,:);
[s1_ss, s2_ss, s3_ss] = size(model);
Tvec_norminal = convert_3d_to_1d(T_3d_nom,s1_ss,s2_ss,s3_ss);
Tvec_abn1 = convert_3d_to_1d(T_3d_abn1,s1_ss,s2_ss,s3_ss);
Tvec_abn2 = convert_3d_to_1d(T_3d_abn2,s1_ss,s2_ss,s3_ss);

load scat_fib_nom_1ghz_down_new.mat;
load scat_fib_nom_1ghz_center_new.mat;
load scat_fib_nom_1ghz_left_new.mat;
load scat_fib_nom_1ghz_right_new.mat;
load scat_fib_nom_1ghz_top_new.mat;
load scat_fib_nom_2ghz_down_new.mat;
load scat_fib_nom_2ghz_center_new.mat;
load scat_fib_nom_2ghz_left_new.mat;
load scat_fib_nom_2ghz_right_new.mat;
load scat_fib_nom_2ghz_top_new.mat;
load scat_fib_nom_3ghz_down_new.mat;
load scat_fib_nom_3ghz_center_new.mat;
load scat_fib_nom_3ghz_left_new.mat;
load scat_fib_nom_3ghz_right_new.mat;
load scat_fib_nom_3ghz_top_new.mat;

%model = model(18:end,:,:);
%[s1_ss, s2_ss, s3_ss] = size(model);
% Read WFs into a cube
wf_cube_1ghz_down = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_1ghz_down_new);
wf_cube_1ghz_center = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_1ghz_center_new);
wf_cube_1ghz_left = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_1ghz_left_new);
wf_cube_1ghz_right = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_1ghz_right_new);
wf_cube_1ghz_top = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_1ghz_top_new);

wf_cube_2ghz_down = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_2ghz_down_new);
wf_cube_2ghz_center = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_2ghz_center_new);
wf_cube_2ghz_left = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_2ghz_left_new);
wf_cube_2ghz_right = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_2ghz_right_new);
wf_cube_2ghz_top = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_2ghz_top_new);

wf_cube_3ghz_down = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_3ghz_down_new);
wf_cube_3ghz_center = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_3ghz_center_new);
wf_cube_3ghz_left = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_3ghz_left_new);
wf_cube_3ghz_right = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_3ghz_right_new);
wf_cube_3ghz_top = read_wf(s1_ss,s2_ss+1,s3_ss,scat_fib_nom_3ghz_top_new);

% eliminate extra data
wf_cube_1ghz_down = wf_cube_1ghz_down(:,1:s2_ss,:);
wf_cube_1ghz_center = wf_cube_1ghz_center(:,1:s2_ss,:);
wf_cube_1ghz_left = wf_cube_1ghz_left(:,1:s2_ss,:);
wf_cube_1ghz_right = wf_cube_1ghz_right(:,1:s2_ss,:);
wf_cube_1ghz_top = wf_cube_1ghz_top(:,1:s2_ss,:);

wf_cube_2ghz_down = wf_cube_2ghz_down(:,1:s2_ss,:);
wf_cube_2ghz_center = wf_cube_2ghz_center(:,1:s2_ss,:);
wf_cube_2ghz_left = wf_cube_2ghz_left(:,1:s2_ss,:);
wf_cube_2ghz_right = wf_cube_2ghz_right(:,1:s2_ss,:);
wf_cube_2ghz_top = wf_cube_2ghz_top(:,1:s2_ss,:);

wf_cube_3ghz_down = wf_cube_3ghz_down(:,1:s2_ss,:);
wf_cube_3ghz_center = wf_cube_3ghz_center(:,1:s2_ss,:);
wf_cube_3ghz_left = wf_cube_3ghz_left(:,1:s2_ss,:);
wf_cube_3ghz_right = wf_cube_3ghz_right(:,1:s2_ss,:);
wf_cube_3ghz_top = wf_cube_3ghz_top(:,1:s2_ss,:);


% find norms of wfs and plot cube results
% convert cube to 1d
wf_vec_1ghz_down = convert_3d_to_1d(wf_cube_1ghz_down,s1_ss,s2_ss,s3_ss);
wf_vec_1ghz_center = convert_3d_to_1d(wf_cube_1ghz_center,s1_ss,s2_ss,s3_ss);
wf_vec_1ghz_left = convert_3d_to_1d(wf_cube_1ghz_left,s1_ss,s2_ss,s3_ss);
wf_vec_1ghz_right = convert_3d_to_1d(wf_cube_1ghz_right,s1_ss,s2_ss,s3_ss);
wf_vec_1ghz_top = convert_3d_to_1d(wf_cube_1ghz_top,s1_ss,s2_ss,s3_ss);

wf_vec_2ghz_down = convert_3d_to_1d(wf_cube_2ghz_down,s1_ss,s2_ss,s3_ss);
wf_vec_2ghz_center = convert_3d_to_1d(wf_cube_2ghz_center,s1_ss,s2_ss,s3_ss);
wf_vec_2ghz_left = convert_3d_to_1d(wf_cube_2ghz_left,s1_ss,s2_ss,s3_ss);
wf_vec_2ghz_right = convert_3d_to_1d(wf_cube_2ghz_right,s1_ss,s2_ss,s3_ss);
wf_vec_2ghz_top = convert_3d_to_1d(wf_cube_2ghz_top,s1_ss,s2_ss,s3_ss);

wf_vec_3ghz_down = convert_3d_to_1d(wf_cube_3ghz_down,s1_ss,s2_ss,s3_ss);
wf_vec_3ghz_center = convert_3d_to_1d(wf_cube_3ghz_center,s1_ss,s2_ss,s3_ss);
wf_vec_3ghz_left = convert_3d_to_1d(wf_cube_3ghz_left,s1_ss,s2_ss,s3_ss);
wf_vec_3ghz_right = convert_3d_to_1d(wf_cube_3ghz_right,s1_ss,s2_ss,s3_ss);
wf_vec_3ghz_top = convert_3d_to_1d(wf_cube_3ghz_top,s1_ss,s2_ss,s3_ss);

% sum and find the norm of each wf
wf_vec_1ghz_norm_down = wf_vec_1ghz_down/sum(wf_vec_1ghz_down);
wf_vec_1ghz_norm_center = wf_vec_1ghz_center/sum(wf_vec_1ghz_center);
wf_vec_1ghz_norm_left   = wf_vec_1ghz_left/sum(wf_vec_1ghz_left);
wf_vec_1ghz_norm_right  = wf_vec_1ghz_right/sum(wf_vec_1ghz_right);
wf_vec_1ghz_norm_top    = wf_vec_1ghz_top/sum(wf_vec_1ghz_top);

wf_vec_2ghz_norm_down = wf_vec_2ghz_down/sum(wf_vec_2ghz_down);
wf_vec_2ghz_norm_center = wf_vec_2ghz_center/sum(wf_vec_2ghz_center);
wf_vec_2ghz_norm_left   = wf_vec_2ghz_left/sum(wf_vec_2ghz_left);
wf_vec_2ghz_norm_right  = wf_vec_2ghz_right/sum(wf_vec_2ghz_right);
wf_vec_2ghz_norm_top    = wf_vec_2ghz_top/sum(wf_vec_2ghz_top);

wf_vec_3ghz_norm_down = wf_vec_3ghz_down/sum(wf_vec_3ghz_down);
wf_vec_3ghz_norm_center = wf_vec_3ghz_center/sum(wf_vec_3ghz_center);
wf_vec_3ghz_norm_left   = wf_vec_3ghz_left/sum(wf_vec_3ghz_left);
wf_vec_3ghz_norm_right  = wf_vec_3ghz_right/sum(wf_vec_3ghz_right);
wf_vec_3ghz_norm_top    = wf_vec_3ghz_top/sum(wf_vec_3ghz_top);


wf_cube_1ghz_norm_down = convert_1d_to_3d(wf_vec_1ghz_norm_down,s1_ss,s2_ss,s3_ss);   % 1ghz
wf_cube_1ghz_norm_center = convert_1d_to_3d(wf_vec_1ghz_norm_center,s1_ss,s2_ss,s3_ss);   % 1ghz
wf_cube_1ghz_norm_left   = convert_1d_to_3d(wf_vec_1ghz_norm_left,s1_ss,s2_ss,s3_ss);     % 1ghz
wf_cube_1ghz_norm_right  = convert_1d_to_3d(wf_vec_1ghz_norm_right,s1_ss,s2_ss,s3_ss);    % 1ghz
wf_cube_1ghz_norm_top    = convert_1d_to_3d(wf_vec_1ghz_norm_top,s1_ss,s2_ss,s3_ss);      % 1ghz

wf_cube_2ghz_norm_down = convert_1d_to_3d(wf_vec_2ghz_norm_down,s1_ss,s2_ss,s3_ss);   % 2ghz
wf_cube_2ghz_norm_center = convert_1d_to_3d(wf_vec_2ghz_norm_center,s1_ss,s2_ss,s3_ss);   % 2ghz
wf_cube_2ghz_norm_left   = convert_1d_to_3d(wf_vec_2ghz_norm_left,s1_ss,s2_ss,s3_ss);     % 2ghz
wf_cube_2ghz_norm_right  = convert_1d_to_3d(wf_vec_2ghz_norm_right,s1_ss,s2_ss,s3_ss);    % 2ghz
wf_cube_2ghz_norm_top    = convert_1d_to_3d(wf_vec_2ghz_norm_top,s1_ss,s2_ss,s3_ss);      % 2ghz

wf_cube_3ghz_norm_down = convert_1d_to_3d(wf_vec_3ghz_norm_down,s1_ss,s2_ss,s3_ss);   % 3ghz
wf_cube_3ghz_norm_center = convert_1d_to_3d(wf_vec_3ghz_norm_center,s1_ss,s2_ss,s3_ss);   % 3ghz
wf_cube_3ghz_norm_left   = convert_1d_to_3d(wf_vec_3ghz_norm_left,s1_ss,s2_ss,s3_ss);     % 3ghz
wf_cube_3ghz_norm_right  = convert_1d_to_3d(wf_vec_3ghz_norm_right,s1_ss,s2_ss,s3_ss);    % 3ghz
wf_cube_3ghz_norm_top    = convert_1d_to_3d(wf_vec_3ghz_norm_top,s1_ss,s2_ss,s3_ss);      % 3ghz

WF = [wf_vec_1ghz_norm_left,wf_vec_2ghz_norm_left,wf_vec_3ghz_norm_left,...
      wf_vec_1ghz_norm_right,wf_vec_2ghz_norm_right,wf_vec_3ghz_norm_right,...
      wf_vec_1ghz_norm_top,wf_vec_2ghz_norm_top,wf_vec_3ghz_norm_top,...
      wf_vec_1ghz_norm_down,wf_vec_2ghz_norm_down,wf_vec_3ghz_norm_down,...
      wf_vec_1ghz_norm_center,wf_vec_2ghz_norm_center,wf_vec_3ghz_norm_center,...
      ];

% Calculate Nominal Brightness Temperature, TB_nominal
%Tvec_norminal = convert_3d_to_1d(T_3d_nom,s1_ss,s2_ss,s3_ss);
TB_norminal = (WF'*Tvec_norminal)';

% Calculate Brightness Temperature with Tumors
%Tvec_abn1 = convert_3d_to_1d(T_3d_abn1,s1_ss,s2_ss,s3_ss);

TB_abn1 = (WF'*Tvec_abn1)';
TB_abn2 = (WF'*Tvec_abn2)';

% Calculate TB_delta
TB_delta_1 = TB_abn1-TB_norminal;
TB_delta_2 = TB_abn2-TB_norminal;

[r_wf, col_wf] = size(WF);
[r_A_norminal, col_A_norminal] = size(A_norminal);
S_norminal = zeros(r_A_norminal, col_wf);
for i = 1:15
    S_norminal(:,i) = cgs(A_norminal',WF(:,i));
end
S_norminal = S_norminal';
% Generate dictionary of tumor position
diameter = 20;
k=1;
i = 1;
ii = 0;
jj = 0;
kk = 0;
for z = 1:20:131
    for x = 21:20:131
        for y = 1:20:141
            tum_col_matrix(i,:) = Gen_Col_Index_Tumor(z,x,y,s1_ss,s2_ss,diameter,k);
            i = i + 1;
            ii = ii + 1;
        end
        jj = jj + 1;
    end
    kk = kk + 1;
end
[r_tum_pos, col_tum_pos] = size(tum_col_matrix);
S_reduce = zeros(col_wf,col_tum_pos);
S_reduce_total = zeros(col_wf,r_tum_pos);
S_reduce_norm_total = zeros(r_tum_pos,1);
TBa_total = zeros(col_wf,r_tum_pos);
TBa_norm_total=zeros(1,r_tum_pos);
for i = 1:r_tum_pos
    for j = 1:col_tum_pos
        S_reduce(:,j) = S_norminal(:,tum_col_matrix(i,j));
    end
    S_reduce = sum(S_reduce,2);
    TBa = S_reduce;
    TBa_norm = norm(TBa);
    TBa_total(:,i) = TBa;
    TBa_norm_total(:,i) = TBa_norm;
end
inner_pro_1 = zeros(r_tum_pos,1);
inner_pro_2 = zeros(r_tum_pos,1);


for i = 1:r_tum_pos
    inner_pro_1(i,1) = vpa(dot(TB_delta_1,TBa_total(:,i)),6)/TBa_norm_total(:,i);
    inner_pro_2(i,1) = vpa(dot(TB_delta_2,TBa_total(:,i)),6)/TBa_norm_total(:,i);
end
[Value_1,Index_1] = max(inner_pro_1);
[Value_2,Index_2] = max(inner_pro_2);

% ROC OF TUMOR AT 1CM
% alpha to threshold test
% alpha = (vpa(dot(TB_delta, TBa_total(:,Index)),6))^2/TBa_norm_total(1,Index);
% plot ROC without tumor, 100 times
% TB_delta_H0 = T_noise;
% TB_delat_H1 = T_noise+c*TBa_hat;
alpha_H0_1 = zeros(100,1);
alpha_H1_1 = zeros(100,1);
T_noise = 0.1*randn(15,100);

TB_delta_H0_1 = T_noise;
TBa_hat_1 = TBa_total(:,Index_1);
TBa_hat_norm_1 =norm(TBa_hat_1);
const_1 = dot(TB_delta_1, TBa_hat_1)/TBa_hat_norm_1^2;
TB_delta_H1_1 = zeros(15,100);
for i = 1:100
    TB_delta_H1_1(:,i) = T_noise(:,i) + const_1*TBa_hat_1;
end

for i = 1:100
    alpha_H0_1(i,1) = (vpa(dot(TB_delta_H0_1(:,i), TBa_hat_1),6))^2/TBa_hat_norm_1^2;
    alpha_H1_1(i,1) = (vpa(dot(TB_delta_H1_1(:,i), TBa_hat_1),6))^2/TBa_hat_norm_1^2;
end    
%t = [3e-8:3e-8:3e-4]; 
t_1 = linspace(min(min(alpha_H0_1),min(alpha_H1_1)),max(max(alpha_H0_1),max(alpha_H1_1)));
[r_t_1, col_t_1] = size(t_1);
N1 = 0;
N2 = 0;
PFA_1 = zeros(1,col_t_1);
PDET_1 = zeros(1,col_t_1);
for i = 1 : col_t_1
    for j = 1:100
        if alpha_H0_1(j,1) > t_1(1,i)
            N2 = N2+1;
        end
    end

    PFA_1(1,i) = N2/100;
    N2 = 0;
   
end
for i = 1:col_t_1
    for j = 1:100
        if alpha_H1_1(j,1) > t_1(1,i)
            N1 = N1 + 1;
        end
    end
    PDET_1(1,i) = N1/100;
    N1 = 0;
end

% ROC of tumor at 2cm 
alpha_H0_2 = zeros(100,1);
alpha_H1_2 = zeros(100,1);
T_noise = 0.1*randn(15,100);

TB_delta_H0_2 = T_noise;
TBa_hat_2 = TBa_total(:,Index_2);
TBa_hat_norm_2 =norm(TBa_hat_2);
const_2 = dot(TB_delta_2, TBa_hat_2)/TBa_hat_norm_2^2;
TB_delta_H1_2 = zeros(15,100);
for i = 1:100
    TB_delta_H1_2(:,i) = T_noise(:,i) + const_2*TBa_hat_2;
end

for i = 1:100
    alpha_H0_2(i,1) = (vpa(dot(TB_delta_H0_2(:,i), TBa_hat_2),6))^2/TBa_hat_norm_2^2;
    alpha_H1_2(i,1) = (vpa(dot(TB_delta_H1_2(:,i), TBa_hat_2),6))^2/TBa_hat_norm_2^2;
end    
t_2 = linspace(min(min(alpha_H0_2),min(alpha_H1_2)),max(max(alpha_H0_2),max(alpha_H1_2)));
[r_t_2, col_t_2] = size(t_2);
N1 = 0;
N2 = 0;
PFA_2 = zeros(1,col_t_2);
PDET_2 = zeros(1,col_t_2);
for i = 1 : col_t_2
    for j = 1:100
        if alpha_H0_2(j,1) > t_2(1,i)
            N2 = N2+1;
        end
    end

    PFA_2(1,i) = N2/100;
    N2 = 0;
   
end
for i = 1:col_t_2
    for j = 1:100
        if alpha_H1_2(j,1) > t_2(1,i)
            N1 = N1 + 1;
        end
    end
    PDET_2(1,i) = N1/100;
    N1 = 0;
end
a = linspace(0,1);
b = linspace(0,1);
figure
plot(PFA_1,PDET_1,'g');hold on;
plot(PFA_2,PDET_2,'r');hold on;
plot(a,b); hold off
xlabel('False Alarm Rate');
ylabel('Detetion Rate');
title('ROC for tumor');
legend('tumor depth is 1 cm','tumor depth is 2 cm','Location','southeast');


   
            


