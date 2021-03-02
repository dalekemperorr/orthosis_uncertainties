%==========================================================================
%[name] two_points_from_vector_xyzijk
%[desc] calculat e position of two arbitrary points located on input axis. 
%       The distance between points is defined by scalling factor 
%[in]   axis_vector -  vector[x y z i j k] describing input axis's orientation
%[in]   scaling_factor - multiplication of distance between two poins
%       defined as length of input's axis <x+y+z> result vector
%[out]  Pa0 - vector [x y z ] describing first point located on input axis
%[out]  Pa1 - vector [x y z ] describing second point located on input axis
%==========================================================================
function [Pa0, Pa1] = two_points_from_vector_xyzijk(axis_vector, scaling_factor)
  
 Pa0 = zeros(3,1);  
 Pa1 = zeros(3,1);  
 
  
  Axis_matrix = eye(4,4);
  Axis_matrix(1,1) = axis_vector(4);
  Axis_matrix(2,1) = axis_vector(5);
  Axis_matrix(3,1) = axis_vector(6);
  Axis_matrix(1,4) = axis_vector(1);
  Axis_matrix(2,4) = axis_vector(2);
  Axis_matrix(3,4) = axis_vector(3);
  
Pa0(1) = axis_vector(1);  
Pa0(2) = axis_vector(2);  
Pa0(3) = axis_vector(3);  
Pa1(1) = axis_vector(1) + axis_vector(4) * scaling_factor;  
Pa1(2) = axis_vector(2) + axis_vector(5) * scaling_factor;  
Pa1(3) = axis_vector(3) + axis_vector(6) * scaling_factor;  

end
