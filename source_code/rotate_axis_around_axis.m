%==========================================================================
%[name] rotate_axis_around_axis
%[desc] rotate one axis around another axis by given angle
%[in]   axis_vector - vector[x y z i j k] describing input axis's orientation 
%[in]   rotation_axis_vector -  vector[x y z i j k] describing rotation axis's orientation
%[in]   angle_rad - angle of rotation in radians
%[out]  axis_vector_rot - vector[x y z i j k] describing rotated axis's orientation
%==========================================================================
function axis_vector_rot = rotate_axis_around_axis(axis_vector, rotation_axis_vector, angle_rad)
  
  [P0, P1] = two_points_from_vector_xyzijk(axis_vector, 5);
  [Pa0, Pa1] = two_points_from_vector_xyzijk(rotation_axis_vector, 5);

  P0_r = rotation_arbitrary_axis(P0, Pa0, Pa1, angle_rad);
  P1_r = rotation_arbitrary_axis(P1, Pa0, Pa1, angle_rad);
  axis_vector_rot = vector_xyzijk_from_two_points(P0_r, P1_r);

end
