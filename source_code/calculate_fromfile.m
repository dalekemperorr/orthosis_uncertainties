%==========================================================================
%[name] calculate_fromfile
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
%[out]  Filename - output filename
%==========================================================================
function [ Delta_Alpha, Delta_Beta, Delta_Gamma, Delta_T, GammaValues_vector] = calculate_fromfile( Axis_X_vector_nom, Axis_Y_vector_nom, Axis_Z_vector_nom, A_kon, B_kon, A_ort, B_ort, S, R, L, Filename)

    %% Orthosis measured angles
    filename = uigetfile('*.csv');
    tic
    InputData = readtable(filename,'Delimiter',';');
    
    steps_nb_table = size(InputData);
    steps_nb = steps_nb_table(1) -2;
    

    %number of combinations of uncertainties under consideration
    Uncert_nb = 512;
    Data = zeros(steps_nb, 33);

    %init vector
    B_rzu_3 = [0, 0, 0]';
    Delta_Alpha = zeros(steps_nb, Uncert_nb);
    Delta_Beta = zeros(steps_nb, Uncert_nb);
    Delta_Gamma = zeros(steps_nb, Uncert_nb);
    Delta_T = zeros(steps_nb, Uncert_nb);

    for i = 1:steps_nb+1
        Alpha_raw = deg2rad(InputData.Alfa(i));
        Beta_raw =  deg2rad(InputData.Beta(i));
        Gamma_raw = deg2rad(InputData.Gamma(i));
        GammaValues_vector(i) = Gamma_raw;            
        %main uncertainities searching loop
        for i_unc = 1:Uncert_nb     
            [Alpha, Beta, Gamma, Axis_X_vector, Axis_Y_vector, Axis_Z_vector] = add_uncertainty(Alpha_raw, Beta_raw, Gamma_raw, Axis_X_vector_nom, Axis_Y_vector_nom, Axis_Z_vector_nom, i_unc -1);

     %       Axis_X_vector = Axis_X_vector_nom;  %debug
     %       Axis_Y_vector = Axis_Y_vector_nom;  %debug
     %       Axis_Z_vector = Axis_Z_vector_nom;  %debug
 
   
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


            %5.43
            Delta_Alpha(i, i_unc) = abs(Alpha_kon - Alpha);
            Delta_Beta(i, i_unc)  = abs(Beta_kon - Beta);
            Delta_Gamma(i, i_unc)  = abs(Gamma_kon - Gamma);
            Delta_A = A_kon_3 - A_ort_3;
            Delta_B = B_kon_3 - B_ort_3;

            %5.44
            Delta_T(i, i_unc)  = (sqrt(sum((A_kon_3-A_ort_3).^2)) + sqrt(sum((B_kon_3-B_ort_3).^2)))/2;
 
            %calculate row index for store partial values in memory
            index = (i-1)*Uncert_nb + i_unc;
            %register all important partial data and save it to .csv file
            Data(index, :) = [Axis_X_vector, Axis_Y_vector, Axis_Z_vector, rad2deg(Alpha_raw), rad2deg(Beta_raw), rad2deg(Gamma_raw),rad2deg(Alpha), rad2deg(Beta), rad2deg(Gamma), rad2deg(Alpha_kon), rad2deg(Beta_kon), rad2deg(Gamma_kon), rad2deg(Alpha_kon - Alpha_raw), rad2deg(Beta_kon - Beta_raw), rad2deg(Gamma_kon - Gamma_raw), sqrt(sum(Delta_A.^2)), sqrt(sum(Delta_B.^2)), (sqrt(sum((A_kon_3-A_ort_3).^2)) + sqrt(sum((B_kon_3-B_ort_3).^2)))/2];
        end        
    end
    
    dlmwrite(Filename, Data(:, 1:33), '-append','delimiter', ';');

end

