%==========================================================================
%[name] angle_correction
%[desc] correct angle to be inside -pi:pi range
%[in]   angle_in - angle in radians
%[out]  angle_out - corrected angle in radians
%==========================================================================
function [angle_out] = angle_correction(angle)

%Correct in loop until angle <-pi, pi>
    while(or(angle < -pi, angle > pi))
        if(angle > pi)
            angle = -2*pi + angle;
        else
            angle = 2*pi + angle;  
        end
    end
    
    angle_out = angle;
end

