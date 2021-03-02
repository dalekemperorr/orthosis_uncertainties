%==========================================================================
%[name] add_uncertainty
%[desc] add uncertainity to position based on iterator
%[in]   V_in - Vector
%[in]   i - iterator, selects proper case for adding uncertainity
%[out]  V_out - Vector
%==========================================================================
function [ VX_out, VY_out, VZ_out ] = add_uncertainty( VX_in, VY_in, VZ_in, i)

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

    %implementation of table 5.29 (combination of uncertainties)
    for j = 6:-1:1
        
        if(not(bitand(mask, i))) 
            %linear
            if(j<=3)
                VX_out(j) = VX_in(j) + U(1,j);
                VY_out(j) = VY_in(j) + U(2,j);
                VZ_out(j) = VZ_in(j) + U(3,j);
                %fprintf('j= %d, +', j)         %debug
            
            %angular - determine sign of uncertainity
            else
                UX(j-3) = U(1,j);
                UY(j-3) = U(2,j);
                UZ(j-3) = U(3,j);                
            end
        else
            %linear
            if(j<=3)
                VX_out(j) = VX_in(j) - U(1,j);
                VY_out(j) = VY_in(j) - U(2,j);
                VZ_out(j) = VZ_in(j) - U(3,j);
                %fprintf('j= %d, -', j)         %debug
            
            %angular - determine sign of uncertainity [deg]
            else
                UX(j-3) = -U(1,j);
                UY(j-3) = -U(2,j);
                UZ(j-3) = -U(3,j);                

            end
          end
        
        %update mask for next vector's item
        mask = bitshift(mask, 1);        
    end
    
    %when all uncertainities signes are determined, calculate actual Axis
    %orientation
    VX_out(4:6) = rotate_single_axis_3D(VX_in, UX);
    VY_out(4:6) = rotate_single_axis_3D(VY_in, UY);
    VZ_out(4:6) = rotate_single_axis_3D(VZ_in, UZ);

end

