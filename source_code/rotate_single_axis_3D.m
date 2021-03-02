%==========================================================================
%[name] rotate_single_axis_3D
%[desc] rotate single axis by 3 given angles
%[in]   axis_vector - vector[x y z i j k] describing input axis's orientation 
%[in]   angle_deg_vector -  vector[a b c] describing rotation axis's in deg
%[out]  axis_vector_rot - vector[x y z i j k] describing rotated axis's orientation
%==========================================================================
function axis_vector_rot = rotate_single_axis_3D(axis_vector, angle_deg_vector)
  
    axis_vector_rot = zeros(3,1);
 
    %Prepare rotation matrixes from uncertainties values
    Rot_X=eye(3);
    Rot_Y=eye(3);
    Rot_Z=eye(3);
    
    Rot_X(2,2) =  cos(deg2rad(angle_deg_vector(1)));
    Rot_X(2,3) = -sin(deg2rad(angle_deg_vector(1)));    
    Rot_X(3,2) =  sin(deg2rad(angle_deg_vector(1)));
    Rot_X(3,3) =  cos(deg2rad(angle_deg_vector(1)));

    Rot_Y(1,1) =  cos(deg2rad(angle_deg_vector(2)));
    Rot_Y(1,3) =  sin(deg2rad(angle_deg_vector(2)));    
    Rot_Y(3,1) = -sin(deg2rad(angle_deg_vector(2)));
    Rot_Y(3,3) =  cos(deg2rad(angle_deg_vector(2)));

    Rot_Z(1,1) =  cos(deg2rad(angle_deg_vector(3)));
    Rot_Z(1,2) = -sin(deg2rad(angle_deg_vector(3)));    
    Rot_Z(2,1) =  sin(deg2rad(angle_deg_vector(3)));
    Rot_Z(2,2) =  cos(deg2rad(angle_deg_vector(3)));
    
    A_rot_x = Rot_X*(axis_vector(4:6)');
    A_rot_y = Rot_Y*(axis_vector(4:6)');
    A_rot_z = Rot_Z*(axis_vector(4:6)');
    
    A_rot_x = A_rot_x./max(A_rot_x);
    A_rot_y = A_rot_y./max(A_rot_y);
    A_rot_z = A_rot_z./max(A_rot_z);
    
    %axis_vector_norm = axis_vector(4:6)'./sqrt(sumsqr(axis_vector(4:6)));
    axis_vector_norm = axis_vector(4:6)'./sqrt(sum(axis_vector(4:6).^2));
    
    A_rot = A_rot_x + A_rot_y + A_rot_z - 2*axis_vector_norm;
    
    axis_vector_rot = A_rot./sqrt(sumsqr(A_rot));
   

end