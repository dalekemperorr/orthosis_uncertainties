%==========================================================================
%[name] add_uncertainty
%[desc] add uncertainity to position based on iterator
%[in]   alpha_in, beta_in, gamma_in - measured angles of orthosis [rad]
%[in]   V_in - Vector xyzijk with Axis orientation
%[in]   i - iterator, selects proper case for adding uncertainity
%[out]  V_out - Vector xyzijk with Axis orientation
%[out]  alpha_out, beta_out, gamma_out - measured angles of orthosis with uncertainties added [rad]
%==========================================================================
function [alpha_out, beta_out, gamma_out, VX_out, VY_out, VZ_out ] = add_uncertainty(alpha_in, beta_in, gamma_in, VX_in, VY_in, VZ_in, i)

    %initialization
    U= zeros(3,6);
     
    mask = 1;
    VX_out = VX_in;
    VY_out = VY_in;
    VZ_out = VZ_in;
 
    %Uncertaintes values [ux mm, uy mm, uz mm, ualfa deg, ubeta deg, ugamma deg]
    U(2,:) = [10.07, 10.55, 10.77, 1.17, 2.58, 1.17];  %Y
    U(3,:) = [14.46, 13.461 16.84, 2.36, 2.33, 2.52];  %Z
    U(1,:) = [12.61, 13.77, 16.83, 2.61, 4.18, 3.22];  %X
   
    %Uncertaintes values for measured angle [uAlphaOrt deg, uBetaOrt, uGammaOrt deg]
    Uort = [U(2,5), U(3,6), U(1,4)];

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
            elseif(j==3)
                gamma_out = gamma_in + deg2rad(Uort(3)); 
            elseif(j==2)
                beta_out  = beta_in  + deg2rad(Uort(2)); 
            elseif(j==1)
                alpha_out = alpha_in + deg2rad(Uort(1)); 
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
            elseif(j==3)
                gamma_out = gamma_in - deg2rad(Uort(3)); 
            elseif(j==2)
                beta_out  = beta_in  - deg2rad(Uort(2)); 
            elseif(j==1)
                alpha_out = alpha_in - deg2rad(Uort(1)); 
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

