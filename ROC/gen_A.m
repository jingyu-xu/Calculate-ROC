function [A] = gen_A(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen)

% tumor_radius = 5; % radius, voxel wise
% tumor_radius_m = tumor_radius;
% tumor_depth = 10; % 
C_qm = 3.27e6;
tau = log(tumor_radius*2/1000/.01)/.002134+50;

% tum_y_cen = 90; %ceil(s2_ss/2);
% tum_x_cen = 20 + tumor_depth;  % skin position + depth + tumor_radius. really is x
% tum_z_cen = floor(s3_ss/2);

bar_norm_prop   = struct('therm_cond',.48,'pcw',2400 ,'qm',700);
bar_gland_prop  = struct('therm_cond',.48,'pcw',2400 ,'qm',700);
bar_fat_prop    = struct('therm_cond',.21,'pcw',800 ,'qm',400);
bar_tumor_prop  = struct('therm_cond',.511,'pcw',48000,'qm', C_qm/tau);  % pcw = 48000, k=.511
bar_air_prop    = struct('therm_cond',0,'pcw',0,'qm',0);
bar_skin_prop   = struct('therm_cond',.343,'pcw',360000,'qm',360); % 1.169e3, 360000, 
bar_muscle_prop = struct('therm_cond', .48,'pcw',3600*1080,'qm',700);
bar_pcw = zeros(s1_ss,s2_ss,s3_ss);
bar_qm  = zeros(s1_ss,s2_ss,s3_ss);
bar_k   = zeros(s1_ss,s2_ss,s3_ss);


W  = zeros(s1_ss,s2_ss,s3_ss);
Qm = zeros(s1_ss,s2_ss,s3_ss);
k  = zeros(s1_ss,s2_ss,s3_ss);
pc = zeros(s1_ss,s2_ss,s3_ss);
model_test = zeros(s1_ss,s2_ss,s3_ss);

for x = 1:s1_ss
    for y = 1:s2_ss
        for z = 1:s3_ss
            id = model(x,y,z);
            model_test(x,y,z) = id;
            if ((( x - tum_x_cen)^2 + (y-tum_y_cen)^2 + (z-tum_z_cen)^2) <= tumor_radius^2) && tumor_on == 1
                prop = bar_tumor_prop;
                model_test(x,y,z) = 5;
            else
                if id == -1
                    prop = bar_air_prop;
                    model_test(x,y,z) = 10;
                elseif id == -2
                    prop = bar_skin_prop;
                elseif id == -4
                    prop = bar_muscle_prop;
%                     prop = bar_norm_prop;
%                     model_test(x,y,z) = -4;
                elseif id == 1.1
                    prop = bar_gland_prop;
%                     prop = bar_norm_prop;
%                     model_test(x,y,z) = 1.1;
                elseif id == 1.2
                    prop = bar_gland_prop;
%                     prop = bar_norm_prop;
%                     model_test(x,y,z) = 1.1;
                elseif id == 1.3
                    prop = bar_gland_prop;
%                     prop = bar_norm_prop;
%                     model_test(x,y,z) = 1.3;
                elseif id == 2
                    prop = bar_norm_prop;
%                     model_test(x,y,z) = 2;
                elseif id == 3.1
                    prop = bar_fat_prop;
%                     prop = bar_norm_prop;
                elseif id == 3.2
                    prop = bar_gland_prop;
%                     prop = bar_norm_prop;
%                     model_test(x,y,z) = 3.2;
                elseif id == 3.3
                    prop = bar_gland_prop;
%                     prop = bar_norm_prop;
%                     model_test(x,y,z) = 3.2;
                end
            end            
            bar_pcw(x,y,z) = prop.pcw;
            bar_k(x,y,z)   = prop.therm_cond;
            bar_qm(x,y,z)  = prop.qm;
        end
    end
end

% figure;contourf(model_test(:,:,floor(s3_ss/2)));

% Tambient = 21; % 27 or 21
% Tart = 37;
A_size = s1_ss*s2_ss*s3_ss;
delta_x = 1/1000;          % convert mm to m, originally .5mm
delta_y = 1/1000;          % with subsamp, 1mm
delta_z = 1/1000;
A_diags = zeros(A_size,7);
A_diags(:,4) = 1;           % don't need to fill main diagonal, its always 1
B = zeros(A_size,1);        % Takes much less time to fill than sparse A matrix
cur_pos = 1;
h = 13.5; %13.5;
use_new = 1;
doing_convec = 1;
for z = 1:s3_ss
    for x = 1:s1_ss %for x = 1:s1_ss
        for y =1:s2_ss %for y = 1:s2_ss
            if model(x,y,z) == -1
                B(cur_pos) = Tambient;
            elseif x > muscle_wall
                B(cur_pos) = Tart;
            elseif model(x,y,z)== -1  || z==1 || x==1 || y==1 || y==s2_ss || z==s3_ss     % air
                if x > skin_start            
                    B(cur_pos) = Tart;  % in body
                else
                    B(cur_pos) = Tambient;
                end
            else
                if use_new == 1
                    % first determine which voxels boarder with air
                    xmin = 0; xplus = 0;
                    ymin = 0; yplus = 0;
                    zmin = 0; zplus = 0; NumB = 0;
                    if doing_convec == 1
                        if (model(x,y,z) == -2 && model(x-1,y,z) == -1)
    %                         keyboard
                            xmin = 1; NumB = NumB + 1;
                        end
                        if model(x,y,z) == -2 && model(x+1,y,z) == -1
                            xplus = 1; NumB = NumB + 1;
                        end
                        if model(x,y,z) == -2 && model(x,y-1,z) == -1
                            ymin = 1; NumB = NumB + 1;
                        end
                        if model(x,y,z) == -2 && model(x,y+1,z) == -1
                            yplus = 1; NumB = NumB + 1;
                        end
                        if model(x,y,z) == -2 && model(x,y,z-1) == -1
                            zmin = 1; NumB = NumB + 1;
                        end
                        if model(x,y,z) == -2 && model(x,y,z+1) == -1
                            zplus = 1; NumB = NumB + 1; 
                        end
                    end

                    wpc = bar_pcw(x,y,z);
                    phi = bar_k(x-1,y,z)/delta_x^2*(~xmin) + bar_k(x+1,y,z)/delta_x^2*(~xplus) + bar_k(x,y-1,z)/delta_y^2*(~ymin) + bar_k(x,y+1,z)/delta_y^2*(~yplus)+ bar_k(x,y,z-1)/delta_z^2*(~zmin) + bar_k(x,y,z+1)/delta_z^2*(~zplus) + wpc;
                    mu = NumB + 1;
                    if xmin == 1
                        tau = h - bar_k(x,y,z)/delta_x;
                        Txmin_coef  = bar_k(x,y,z)/delta_x/tau/mu;
                    else
                        Txmin_coef  = -bar_k(x-1,y,z)/delta_x^2/phi/mu;
                    end
                    if xplus == 1
                        tau = h - bar_k(x,y,z)/delta_x;
                        Txplus_coef  = bar_k(x,y,z)/delta_x/tau/mu;
                    else
                        Txplus_coef  = -bar_k(x+1,y,z)/delta_x^2/phi/mu;
                    end
                    if model(x,y,z) == -2 && model(x,y-1,z) == -1
                        tau = h - bar_k(x,y,z)/delta_y;
                        Tymin_coef  =  bar_k(x,y,z)/delta_x/tau/mu;
                    else
                        Tymin_coef  = -bar_k(x,y-1,z)/delta_y^2/phi/mu;
                    end
                    if model(x,y,z) == -2 && model(x,y+1,z) == -1
                        tau = h - bar_k(x,y,z)/delta_y;
                        Typlus_coef  = bar_k(x,y,z)/delta_x/tau/mu;
                    else
                        Typlus_coef  = -bar_k(x,y+1,z)/delta_y^2/phi/mu;
                    end
                    if model(x,y,z) == -2 && model(x,y,z-1) == -1
                        tau = h - bar_k(x,y,z)/delta_z;
                        Tzmin_coef  = bar_k(x,y,z)/delta_x/tau/mu;
                    else
                        Tzmin_coef  = -bar_k(x,y,z-1)/delta_z^2/phi/mu;
                    end
                    if model(x,y,z) == -2 && model(x,y,z+1) == -1
                        tau = h - bar_k(x,y,z)/delta_z;
                        Tzplus_coef  =  bar_k(x,y,z)/delta_x/tau/mu;
                    else
                        Tzplus_coef  = -bar_k(x,y,z+1)/delta_z^2/phi/mu;
                    end
                    
%                     Txmin_coef   = -bar_k(x-1,y,z)/delta_x^2/phi;
%                     Txplus_coef  = -bar_k(x+1,y,z)/delta_x^2/phi;
%                     Tymin_coef   = -bar_k(x,y-1,z)/delta_y^2/phi;
%                     Typlus_coef  = -bar_k(x,y+1,z)/delta_y^2/phi;
%                     Tzmin_coef   = -bar_k(x,y,z-1)/delta_z^2/phi;
%                     Tzplus_coef  = -bar_k(x,y,z+1)/delta_z^2/phi;

                    A_diags(cur_pos-s2_ss,2)       = Txmin_coef;
                    A_diags(cur_pos+s2_ss,6)       = Txplus_coef;
                    A_diags(cur_pos-1,3)           = Tymin_coef;
                    A_diags(cur_pos+1,5)           = Typlus_coef;
                    A_diags(cur_pos-s1_ss*s2_ss,1) = Tzmin_coef; 
                    A_diags(cur_pos+s1_ss*s2_ss,7) = Tzplus_coef;
%                     B(cur_pos) = (bar_qm(x,y,z) + wpc*Tart + NumB*h/(h-bar_k(x,y,z)/delta_x)*Tambient)/phi;
                    
                    B(cur_pos) = (bar_qm(x,y,z) + wpc*Tart)/phi/mu + h/tau*Tambient*NumB/mu*doing_convec;
                else
                    convec_used = 0;  % reset varrs
                    x_convec_on = 0;   %#ok<*NASGU>
                    y_convec_on = 0;
                    z_convec_on = 0;
                    if model(x,y,z) == -50
                        if (model(x-1,y,z) == -1 || model(x+1,y,z) == -1)
                            tau = h+bar_k(x,y,z)/delta_x;
                            A_diags(cur_pos-s2_ss,2) = bar_k(x,y,z)/delta_x/tau; %-bar_k(x,y,z)/h;
                            A_diags(cur_pos+s2_ss,6) = bar_k(x,y,z)/delta_x/tau; %-bar_k(x,y,z)/h;
                            convec_used = 1;
                            B(cur_pos) = Tambient;
%                             B(cur_pos) = h*Tambient/tau;
                            x_convec_on = 1;
                        end
                        if (model(x,y-1,z) == -1 || model(x,y+1,z) == -1)
                            tau = h+bar_k(x,y,z)/delta_y;
                            A_diags(cur_pos-1,3) = bar_k(x,y,z)/delta_y/tau; %-bar_k(x,y,z)/h;
                            A_diags(cur_pos+1,5) = bar_k(x,y,z)/delta_y/tau; %-bar_k(x,y,z)/h;
                            convec_used = 1;
%                             B(cur_pos) = h*Tambient/tau;
                            B(cur_pos) = Tambient;
                            y_convec_on = 1;
                        end
                        if (model(x,y,z-1) == -1 || model(x,y,z+1) == -1)
                            tau = h+bar_k(x,y,z)/delta_z;
                            A_diags(cur_pos-s1_ss*s2_ss,1) = bar_k(x,y,z)/delta_z/tau; %-bar_k(x,y,z)/h;
                            A_diags(cur_pos+s1_ss*s2_ss,7) = bar_k(x,y,z)/delta_z/tau; %-bar_k(x,y,z)/h;
                            convec_used = 1;
                            B(cur_pos) = Tambient;
%                             B(cur_pos) = h*Tambient/tau;
                            z_convec_on = 1;
                        end
                    end
                    if convec_used == 0
%                         keyboard
                        wpc = bar_pcw(x,y,z);
                        phi = bar_k(x-1,y,z)/delta_x^2 + bar_k(x+1,y,z)/delta_x^2 + bar_k(x,y-1,z)/delta_y^2 + bar_k(x,y+1,z)/delta_y^2+ bar_k(x,y,z-1)/delta_z^2 + bar_k(x,y,z+1)/delta_z^2;

%                         bar_k_sum_x = (bar_k(x,y,z)+bar_k(x,y,z))/delta_x^2;
%                         bar_k_sum_y = (bar_k(x,y,z)+bar_k(x,y,z))/delta_y^2;
%                         bar_k_sum_z = (bar_k(x,y,z)+bar_k(x,y,z))/delta_z^2; 2880000
                        bar_k_sum_x = (bar_k(x-1,y,z)+bar_k(x+1,y,z))/delta_x^2;
                        bar_k_sum_y = (bar_k(x,y-1,z)+bar_k(x,y+1,z))/delta_y^2;
                        bar_k_sum_z = (bar_k(x,y,z-1)+bar_k(x,y,z+1))/delta_z^2;
                        k_sum = bar_k_sum_x + bar_k_sum_y + bar_k_sum_z;
%                         beta = 1/k_sum;
                        beta = 1/phi;
                        
                        Txmin_coef   = -bar_k(x-1,y,z)/delta_x^2/phi;
                        Txplus_coef  = -bar_k(x+1,y,z)/delta_x^2/phi;
                        Tymin_coef   = -bar_k(x,y-1,z)/delta_y^2/phi;
                        Typlus_coef  = -bar_k(x,y+1,z)/delta_y^2/phi;
                        Tzmin_coef   = -bar_k(x,y,z-1)/delta_z^2/phi;
                        Tzplus_coef  = -bar_k(x,y,z+1)/delta_z^2/phi;
                        
                        gamma = 1/(1+wpc*beta);
                        test1 = 1/(phi + wpc);
                        test2 = gamma*beta;

                        A_diags(cur_pos-s1_ss*s2_ss,1) = -bar_k(x,y,z-1)/delta_z^2*beta*gamma; 
                        A_diags(cur_pos-s2_ss,2)       = -bar_k(x-1,y,z)/delta_x^2*beta*gamma;
                        A_diags(cur_pos-1,3)           = -bar_k(x,y-1,z)/delta_y^2*beta*gamma;
                        A_diags(cur_pos+1,5)           = -bar_k(x,y+1,z)/delta_y^2*beta*gamma;
                        A_diags(cur_pos+s2_ss,6)       = -bar_k(x+1,y,z)/delta_x^2*beta*gamma;
                        A_diags(cur_pos+s1_ss*s2_ss,7) = -bar_k(x,y,z+1)/delta_z^2*beta*gamma;
                        B(cur_pos) = (bar_qm(x,y,z) + wpc*Tart)*beta*gamma;
%                         B(cur_pos) = (bar_qm(x,y,z) + wpc*Tart)/phi;
                        
%                         A_diags(cur_pos-s1_ss*s2_ss,1) = Tzmin_coef; %-bar_k(x,y,z-1)*beta/delta_z^2*gamma; 
%                         A_diags(cur_pos-s2_ss,2)       = Txmin_coef; %-bar_k(x-1,y,z)*beta/delta_x^2*gamma;
%                         A_diags(cur_pos-1,3)           = Tymin_coef; %-bar_k(x,y-1,z)*beta/delta_y^2*gamma;
%                         A_diags(cur_pos+1,5)           = Typlus_coef; %-bar_k(x,y+1,z)*beta/delta_y^2*gamma;
%                         A_diags(cur_pos+s2_ss,6)       = Txplus_coef; %-bar_k(x+1,y,z)*beta/delta_x^2*gamma;
%                         A_diags(cur_pos+s1_ss*s2_ss,7) = Tzplus_coef; %-bar_k(x,y,z+1)*beta/delta_z^2*gamma;
%                         B(cur_pos) = (bar_qm(x,y,z) + wpc*Tart)/phi;
                    end
                end    
            end  
            cur_pos = cur_pos + 1;
        end
    end
end

d = [-s1_ss*s2_ss,-s2_ss,-1,0,1,s2_ss,s1_ss*s2_ss];
A = spdiags(A_diags,d,A_size,A_size);




