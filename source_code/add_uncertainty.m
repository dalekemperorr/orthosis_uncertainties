%==========================================================================
%[name] add_uncertainty
%[desc] add uncertainity to position based on iterator
%[in]   V_in - Vector
%[in]   i - iterator, selects proper case for adding uncertainity
%[out]  V_out - Vector
%==========================================================================
function [alpha, beta, gamma, VX_out, VY_out, VZ_out ] = add_uncertainty(alpha, beta, gamma, VX_in, VY_in, VZ_in, i)

    %initialization
    U= zeros(3,6);
     
    mask = 1;
    VX_out = VX_in;
    VY_out = VY_in;
    VZ_out = VZ_in;
 
    %Uncertaintes values [ux mm, uy mm, uz mm, ualfa deg, ubeta deg, ugamma deg]
    U(1,:) = [17.61, 17.98, 23.06, 4.61, 5.76, 4.74];  %X
    U(2,:) = [15.23, 15.55, 15.77, 4.17, 2.58, 2.17];  %Y
    U(3,:) = [20.48, 18.46, 21.84, 5.36, 3.91, 2.52];  %Z
    
    %Uncertaintes values for measured angle [uAlphaOrt deg, uBetaOrt, uGammaOrt deg]
    Uort = [5.0, 5.0, 9.0];

    %implementation of table 5.29 (combination of uncertainties)
    for j = 9:-1:1
        
        if(not(bitand(mask, i))) 
           %angular - determine sign of uncertainity
            if(j>6)
                UX(j-6) = U(1,j-3);
                UY(j-6) = U(2,j-3);
                UZ(j-6) = U(3,j-3);                
            %linear
            elseif(and(j<=6, j>3))
                VX_out(j-3) = VX_in(j-3) + U(1,j-3);
                VY_out(j-3) = VY_in(j-3) + U(2,j-3);
                VZ_out(j-3) = VZ_in(j-3) + U(3,j-3);
                %fprintf('j= %d, +', j)         %debug
            %orthosis measured angle
            else
                alpha = alpha + deg2rad(Uort(1)); 
                beta  = beta  + deg2rad(Uort(2)); 
                gamma = gamma + deg2rad(Uort(3)); 
             end
                
            
         else
            %angular - determine sign of uncertainity
            if(j>6)
                UX(j-6) = -U(1,j-3);
                UY(j-6) = -U(2,j-3);
                UZ(j-6) = -U(3,j-3);                
            %linear
            elseif(and(j<=6, j>3))
                VX_out(j-3) = VX_in(j-3) - U(1,j-3);
                VY_out(j-3) = VY_in(j-3) - U(2,j-3);
                VZ_out(j-3) = VZ_in(j-3) - U(3,j-3);
                %fprintf('j= %d, +', j)         %debug
            %orthosis measured angle
            else
                alpha = alpha - deg2rad(Uort(1)); 
                beta  = beta  - deg2rad(Uort(2)); 
                gamma = gamma - deg2rad(Uort(3)); 
            end
         end
        
        %update mask for next vector's item
        mask = bitshift(mask, 1, 'uint16');        
    end
    
    %when all uncertainities signes are determined, calculate actual Axis
    %orientation
    VX_out(4:6) = rotate_single_axis_3D(VX_in, UX);
    VY_out(4:6) = rotate_single_axis_3D(VY_in, UY);
    VZ_out(4:6) = rotate_single_axis_3D(VZ_in, UZ);
end

