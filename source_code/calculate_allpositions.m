%==========================================================================
%[name] calculate_allpositions
%[desc] calculates errors for all combinations of uncertainties
%[in]   Axis_X_vector_nom - Vector xyzijk describing position and
%orientation of Orthosis X axis
%[in]   Axis_Y_vector_nom - Vector xyzijk describing position and
%orientation of Orthosis Y axis
%[in]   Axis_Z_vector_nom - Vector xyzijk describing position and
%orientation of Orthosis Z axis
%[in]   A_kon - position of point A placed on limb wrt global RF
%[in]   B_kon - position of point B placed on limb wrt global RF
%[in]   A_ort - position of point A placed on orthosis wrt global RF
%[in]   B_ort - position of point B placed on orthosis wrt global RF
%[out]  Delta_Alpha - 4D matrix with uncerainties for each alpha, beta,
%gamma combination
%[out]  Delta_Betha - 4D matrix with uncerainties for each alpha, beta,
%gamma combination
%[out]  Delta_Gamma - 4D matrix with uncerainties for each alpha, beta,
%gamma combination
%[out]  Delta_T - 4D matrix with uncerainties for each alpha, beta,
%gamma combination
%==========================================================================
function [ Delta_Alpha, Delta_Beta, Delta_Gamma, Delta_T, GammaValues_vector ] = calculate_allpositions( Axis_X_vector_nom, Axis_Y_vector_nom, Axis_Z_vector_nom, A_kon, B_kon, A_ort, B_ort, S, R, L, Filename )

    tic
 %   cluster = parcluster;
 %   cluster.NumWorkers = 8;
    parallerCalc = parpool(4)
    %% Orthosis measured angles
    %Alpha [deg]
    %Alpha_step = 10;
    %Alpha_start = -30;
    %Alpha_end = 27;
    %determine number of steps
    Alpha_step = 5;
    Alpha_start = -20;
    Alpha_end = 20;
     Alpha_steps_nb = floor((Alpha_end - Alpha_start)/Alpha_step) +1;
    if(Alpha_end > Alpha_start+ Alpha_step*Alpha_steps_nb) 
        Alpha_steps_nb = Alpha_steps_nb+1;
    end
    if(not(mod(Alpha_start, Alpha_step)))
        Alpha_steps_nb = Alpha_steps_nb-1;
    end
    %Beta [deg]
    %Beta_step = 10;
    %Beta_start = -40;
    %Beta_end = 38;
    Beta_step = 5;
    Beta_start = -35;
    Beta_end = 35;
    %determine number of steps
    Beta_steps_nb = floor((Beta_end - Beta_start)/Beta_step) +1;
    if(Beta_end > Beta_start+ Beta_step*Beta_steps_nb) 
        Beta_steps_nb = Beta_steps_nb+1;
    end
    if(not(mod(Beta_start, Beta_step)))
        Beta_steps_nb = Beta_steps_nb-1;
    end
    %Gamma [deg]
   % Gamma_step = 10;
    %Gamma_start = -10;
    %Gamma_end = 120;
    Gamma_step = 5;
    Gamma_start = -5;
    Gamma_end = 85;
    %determine number of steps
    Gamma_steps_nb = floor((Gamma_end - Gamma_start)/Gamma_step) +1;
    if(Gamma_end > floor(Gamma_start+ Gamma_step*Gamma_steps_nb)) 
        Gamma_steps_nb = Gamma_steps_nb+1;
    end
    if(not(mod(Gamma_start, Gamma_step)))
        Gamma_steps_nb = Gamma_steps_nb-1;
    end
    
    GammaValues_vector = zeros(1,Gamma_steps_nb);
 
    %number of combinations of uncertainties under consideration
    Uncert_nb = 512;
    
    Data = zeros(Alpha_steps_nb*Beta_steps_nb*Gamma_steps_nb*Uncert_nb, 33);
    Data_temp = zeros(Uncert_nb, 33);

    %init vector
    B_rzu_3 = [0, 0, 0]';
    Delta_Alpha = zeros(Alpha_steps_nb, Beta_steps_nb, Gamma_steps_nb, Uncert_nb);
    Delta_Beta = zeros(Alpha_steps_nb, Beta_steps_nb, Gamma_steps_nb, Uncert_nb);
    Delta_Gamma = zeros(Alpha_steps_nb, Beta_steps_nb, Gamma_steps_nb, Uncert_nb);
    Delta_T = zeros(Alpha_steps_nb, Beta_steps_nb, Gamma_steps_nb, Uncert_nb);

    %Now real calculations: starting with 4 nested loops- 3 for measured angles
    %sweep and 4th for uncertainties combination
    for i_alpha = 1:Alpha_steps_nb
        if (i_alpha ==1)
            Alpha_raw = Alpha_start;
            Alpha_raw = deg2rad(Alpha_start);
        elseif (i_alpha == Alpha_steps_nb)
            Alpha_raw = Alpha_end;
            Alpha_raw = deg2rad(Alpha_end); 
        else
           Alpha_raw = Alpha_start - mod(Alpha_start, Alpha_step) + (i_alpha -1)*Alpha_step;
           Alpha_raw = deg2rad(Alpha_start - mod(Alpha_start, Alpha_step) + (i_alpha -1)*Alpha_step);    
        end

    for i_beta = 1:Beta_steps_nb
         if (i_beta ==1)
            Beta_raw = Beta_start;
            Beta_raw = deg2rad(Beta_start);
        elseif (i_beta == Beta_steps_nb)
            Beta_raw = Beta_end;
            Beta_raw = deg2rad(Beta_end); 
        else
           Beta_raw = Beta_start - mod(Beta_start, Beta_step) + (i_beta -1)*Beta_step;
           Beta_raw = deg2rad(Beta_start - mod(Beta_start, Beta_step) + (i_beta -1)*Beta_step);    
         end

    for i_gamma = 1:Gamma_steps_nb
        if (i_gamma ==1)
            Gamma_raw = Gamma_start;
            Gamma_raw = deg2rad(Gamma_start);
        elseif (i_gamma == Gamma_steps_nb)
            Gamma_raw = Gamma_end;
            Gamma_raw = deg2rad(Gamma_end); 
        else
           Gamma_raw = Gamma_start - mod(Gamma_start, Gamma_step) + (i_gamma -1)*Gamma_step;
           Gamma_raw = deg2rad(Gamma_start - mod(Gamma_start, Gamma_step) + (i_gamma -1)*Gamma_step);    
        end
        if(and(i_alpha==1, i_beta == 1))
            GammaValues_vector(i_gamma) = rad2deg(Gamma_raw);
        end
       %main uncertainities searching loop
        parfor i_unc = 1:Uncert_nb 

 %           Data_temp = zeros(Uncert_nb, 30);
            [Alpha, Beta, Gamma, Axis_X_vector, Axis_Y_vector, Axis_Z_vector] = add_uncertainty(Alpha_raw, Beta_raw, Gamma_raw, Axis_X_vector_nom, Axis_Y_vector_nom, Axis_Z_vector_nom, i_unc -1);
            
            %Axis_X_vector = Axis_X_vector_nom;  %debug
            %Axis_Y_vector = Axis_Y_vector_nom;  %debug
            %Axis_Z_vector = Axis_Z_vector_nom;  %debug

            %5.23
            A_ort_3 = rotate_point_around_3tilted_axis(A_ort, Axis_X_vector, Axis_Y_vector, Axis_Z_vector, Alpha, Beta, Gamma);
            %5.24
            B_ort_3 = rotate_point_around_3tilted_axis(B_ort, Axis_X_vector, Axis_Y_vector, Axis_Z_vector, Alpha, Beta, Gamma);

            SA_ort_3 = A_ort_3 - S;

            %5.25
            SA_kon_3 = SA_ort_3/sqrt(sum(SA_ort_3.^2))*R;
            A_kon_3 = SA_kon_3;

            %5.32
            L1_t = SA_kon_3(1)*B_ort_3(1)+SA_kon_3(2)*B_ort_3(2)+SA_kon_3(3)*B_ort_3(3);
            L2_t = SA_kon_3(1)*A_kon_3(1)+SA_kon_3(2)*A_kon_3(2)+SA_kon_3(3)*A_kon_3(3);
            M_t = SA_kon_3(1)*SA_kon_3(1)+SA_kon_3(2)*SA_kon_3(2)+SA_kon_3(3)*SA_kon_3(3);
            t = -(L1_t-L2_t)/M_t;

            %5.33
            B_rzu_3 = B_ort_3 + SA_kon_3*t;


            %5.34
            A_kon_3_B_kon_3 = L*(B_rzu_3 - A_kon_3)/sqrt(sum((B_rzu_3 - A_kon_3).^2));
            B_kon_3 = A_kon_3_B_kon_3 + A_kon_3;


            %5.35
            k_3 = (A_kon_3 - B_kon_3)/sqrt(sum((A_kon_3 - B_kon_3).^2));

            %5.36
            SA_kon_3 = A_kon_3 - S;
            i_3 = SA_kon_3/sqrt(sum(SA_kon_3.^2));

            %5.37
            j_3 = cross(k_3, i_3);

            %5.38
            R_kon_3 = [i_3, j_3, k_3];

            %5.40
            Beta_kon_v1 = angle_correction(asin(R_kon_3(2,1)));
            Beta_kon_v2 = angle_correction(pi - asin(R_kon_3(2,1)));

            %5.41
            Gamma_kon_v1 = atan2(-R_kon_3(2,3)/cos(Beta_kon_v1), R_kon_3(2,2)/cos(Beta_kon_v1));
            Gamma_kon_v2 = atan2(-R_kon_3(2,3)/cos(Beta_kon_v2), R_kon_3(2,2)/cos(Beta_kon_v2));

            %5.42
            Alpha_kon_v1 = atan2(-R_kon_3(3,1)/cos(Beta_kon_v1), R_kon_3(1,1)/cos(Beta_kon_v1));
            Alpha_kon_v2 = atan2(-R_kon_3(3,1)/cos(Beta_kon_v2), R_kon_3(1,1)/cos(Beta_kon_v2));


            %%choose proper angle set
            error_v1 = abs((Alpha - Alpha_kon_v1)) + abs((Beta - Beta_kon_v1)) + abs((Gamma - Gamma_kon_v1));
            error_v2 = abs((Alpha - Alpha_kon_v2)) + abs((Beta - Beta_kon_v2)) + abs((Gamma - Gamma_kon_v2));
            if(error_v1 < error_v2)
                Alpha_kon = Alpha_kon_v1;
                Beta_kon  = Beta_kon_v1;
                Gamma_kon = Gamma_kon_v1;    
            else
                Alpha_kon = Alpha_kon_v2;
                Beta_kon  = Beta_kon_v2;
                Gamma_kon = Gamma_kon_v2;      
            end
            
            %force limits on angles
            [Alpha_kon, Beta_kon, Gamma_kon] = force_limits(Alpha_kon, Beta_kon, Gamma_kon);

 
            %5.43
            Delta_Alpha(i_alpha, i_beta, i_gamma, i_unc) = abs(Alpha_kon - Alpha);
            Delta_Beta(i_alpha, i_beta, i_gamma, i_unc)  = abs(Beta_kon - Beta);
            Delta_Gamma(i_alpha, i_beta, i_gamma, i_unc)  = abs(Gamma_kon - Gamma);
            Delta_A = A_kon_3 - A_ort_3;
            Delta_B = B_kon_3 - B_ort_3;

            %5.44
            Delta_T(i_alpha, i_beta, i_gamma, i_unc)  = (sqrt(sum((A_kon_3-A_ort_3).^2)) + sqrt(sum((B_kon_3-B_ort_3).^2)))/2;

            %register all important partial data for further save it to .csv file
            Data_temp(i_unc, :) = [Axis_X_vector, Axis_Y_vector, Axis_Z_vector, rad2deg(Alpha_raw), rad2deg(Beta_raw), rad2deg(Gamma_raw),rad2deg(Alpha), rad2deg(Beta), rad2deg(Gamma), rad2deg(Alpha_kon), rad2deg(Beta_kon), rad2deg(Gamma_kon), rad2deg(Alpha_kon - Alpha_raw), rad2deg(Beta_kon - Beta_raw), rad2deg(Gamma_kon - Gamma_raw), sqrt(sum(Delta_A.^2)), sqrt(sum(Delta_B.^2)), (sqrt(sum((A_kon_3-A_ort_3).^2)) + sqrt(sum((B_kon_3-B_ort_3).^2)))/2];
        end
         
    %progress display
    index = (i_alpha-1)*Beta_steps_nb*Gamma_steps_nb*Uncert_nb + (i_beta-1)*Gamma_steps_nb*Uncert_nb + (i_gamma-1)*Uncert_nb + 1;
    steps_nb = index + Uncert_nb-1;
    Progress_string = sprintf('Calculated %d / %d steps(%.2f %%)', steps_nb, Alpha_steps_nb*Beta_steps_nb*Gamma_steps_nb*Uncert_nb, steps_nb/(Alpha_steps_nb*Beta_steps_nb*Gamma_steps_nb*Uncert_nb)*100);
    clc
    disp(Progress_string)
    %copy from temporary variable for further usage (save to file)
    Data(index:index+Uncert_nb-1, :) = Data_temp;
  end
 
    end
    end
    delete(gcp('nocreate'))
    dlmwrite(Filename, Data(:, 1:33), '-append','delimiter', ';'); 

end

