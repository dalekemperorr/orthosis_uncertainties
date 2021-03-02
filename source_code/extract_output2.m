%==========================================================================
%[name] extract_output
%[desc] extracts multidimentional uncertainties into digestable tables form
%[in]  Delta_Alpha - 2D/4D matrix with uncerainties for each alpha, beta,
%gamma combination
%[in]  Delta_Betha - 2D/4D matrix with uncerainties for each alpha, beta,
%gamma combination
%[in]  Delta_Gamma - 2D/4D matrix with uncerainties for each alpha, beta,
%gamma combination
%[in]  Delta_T - 2D/4D matrix with uncerainties for each alpha, beta,
%gamma combination
function [T] = extract_output2( Delta_Alpha, Delta_Beta, Delta_Gamma, Delta_T, GammaValues_vector )

Table_size = size(Delta_Alpha);

Table_Dimention = size(Table_size);

if (Table_Dimention(2) == 2)
    %handling for angles read from file       
    m = Table_size(1);
    T = zeros (m, 4);

    for(i=1:m)
        T(i, 1)= rad2deg(max(Delta_Alpha(i,:)));
        T(i, 2)= rad2deg(max(Delta_Beta(i,:)));
        T(i, 3)= rad2deg(max(Delta_Gamma(i,:)));
        T(i, 4)= max(Delta_T(i,:));

    end
    
elseif (Table_Dimention(2) == 4)
    %handling for all possible angles  
    m = Table_size(1);
    n = Table_size(2);
    T = zeros (2*m, 4*n);

    for(i=1:m)
        for(j = 1:n)
            Gmax = max(squeeze(Delta_Alpha(i,j,:,:))');
            [Fi, I]= max(Gmax);
            T(2*i-1, 4*j-3) = rad2deg(Fi);
            %index of gamma, need full value
            T(2*i-1, 4*j-2) = GammaValues_vector(I);

            Gmax = max(squeeze(Delta_Beta(i,j,:,:))');
            [Fi, I]= max(Gmax);
            T(2*i-1, 4*j-1) = rad2deg(Fi);
            %index of gamma, need full value
            T(2*i-1, 4*j) = GammaValues_vector(I);

            Gmax = max(squeeze(Delta_Gamma(i,j,:,:))');
            [Fi, I]= max(Gmax);
            T(2*i, 4*j-3) = rad2deg(Fi);
            %index of gamma, need full value
            T(2*i, 4*j-2) = GammaValues_vector(I);

            Gmax = max(squeeze(Delta_T(i,j,:,:))');
            [T(2*i, 4*j-1), I]= max(Gmax);
            %index of gamma, need full value
            T(2*i, 4*j) = GammaValues_vector(I);

        end
    end

   

end
end

