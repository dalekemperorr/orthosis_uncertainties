%==========================================================================
%[name] force_limits
%[desc] forces angles to fit inside arbitrary limits
%[in]   alpha_in - angle in radians
%[in]   beta_in - angle in radians
%[in]   gamma_in - angle in radians
%[out]  alpha_out - corrected angle in radians
%[out]  beta_out - corrected angle in radians
%[out]  gamma_out - corrected angle in radians
%==========================================================================
function [ alpha_out, beta_out, gamma_out ] = force_limits( alpha_in, beta_in, gamma_in )

%Arbitrary limits [limit_min, limit_max] for each angle in degrees
alpha_limits_deg = [ -5, 40];
beta_limits_deg = [ -5, 40];
gamma_limits_deg = [ -5, 90];

%should be calculated once and static for better performance
alpha_limits_rad = deg2rad(alpha_limits_deg);
beta_limits_rad = deg2rad(beta_limits_deg);
gamma_limits_rad = deg2rad(gamma_limits_deg);

%check limits of alpha
    if(alpha_in < alpha_limits_rad(1))
        alpha_out = alpha_limits_rad(1);
    elseif (alpha_in > alpha_limits_rad(2))
        alpha_out = alpha_limits_rad(2);
    else
        alpha_out = alpha_in;
    end
    
%check limits of beta
   if(beta_in < beta_limits_rad(1))
        beta_out = beta_limits_rad(1);
    elseif (beta_in > beta_limits_rad(2))
        beta_out = beta_limits_rad(2);
    else
        beta_out = beta_in;
    end
        
%check limits of gamma
   if(gamma_in < gamma_limits_rad(1))
        gamma_out = gamma_limits_rad(1);
    elseif (gamma_in > gamma_limits_rad(2))
        gamma_out = gamma_limits_rad(2);
    else
        gamma_out = gamma_in;
    end

end

