%==========================================================================
%[name] rotate_point_around_axis
%[desc] rotate point around any axis by given angle 
%[in]   P - veector [x y z] describing input point's position 
%[in]   rotation_axis_vector -  vector[x y z i j k] describing rotation axis's orientation
%[in]   angle_rad - angle of rotation in radians
%[out]  P_rot - vector[x y z ] describing rotated point's orientation
%==========================================================================
function P_rot = rotate_point_around_axis(P, rotation_axis_vector, angle_rad)

 [Pa0, Pa1] = two_points_from_vector_xyzijk(rotation_axis_vector, 5);
  
 P_rot = rotation_arbitrary_axis(P, Pa0, Pa1, angle_rad);
 
end