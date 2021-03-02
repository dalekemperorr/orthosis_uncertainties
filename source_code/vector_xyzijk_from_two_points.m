%==========================================================================
%[name] vector_xyzijk_from_two_points
%[desc] Generate  vector [x y z i j k] that fully describes axis build on
%       two input poins
%[in]  Pa0 - vector [x y z ] describing first point position
%[in]  Pa1 - vector [x y z ] describing second point position
%[in]  axis_vector -  vector[x y z i j k] describing axis buils on input
%points
%==========================================================================
function axis_vector = vector_xyzijk_from_two_points(P0, P1)
  
 L = sqrt((P1(1) - P0(1))^2 + (P1(2) - P0(2))^2 + (P1(3) - P0(3))^2);
axis_vector = zeros(1,6);

axis_vector(1) = P0(1); 
axis_vector(2) = P0(2); 
axis_vector(3) = P0(3); 
axis_vector(4) = (P1(1) - P0(1))/L; 
axis_vector(5) = (P1(2) - P0(2))/L; 
axis_vector(6) = (P1(3) - P0(3))/L; 
 
 
end
