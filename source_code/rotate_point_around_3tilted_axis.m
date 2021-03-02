%==========================================================================
%[name] rotate_point_around_3tilted_axis
%[desc] rotate point around 3 nonortogonal axis by 3 angles
%[in]   P - vector[x y z] describing position of input point 
%[in]   Axis_X_vector - vector[x y z i j k] describing input axis X's  orientation 
%[in]   Axis_Y_vector - vector[x y z i j k] describing input axis Y's  orientation 
%[in]   Axis_Z_vector - vector[x y z i j k] describing input axis Z's  orientation 
%[in]   Alpha -  angle of rotation around Y axis in radians
%[in]   Beta  -  angle of rotation around Z axis in radians
%[in]   Gamma -  angle of rotation around X axis in radians
%[out]  P_3 - vector[x y z] describing position of triple rotated point
%==========================================================================
function [P_3] = rotate_point_around_3tilted_axis(P, Axis_X_vector, Axis_Y_vector, Axis_Z_vector, Alpha, Beta, Gamma )

%%Step 1
P_1 = rotate_point_around_axis(P, Axis_Y_vector, Alpha);
Axis_X_1_vector = rotate_axis_around_axis(Axis_X_vector,Axis_Y_vector,Alpha);
Axis_Z_1_vector = rotate_axis_around_axis(Axis_Z_vector,Axis_Y_vector,Alpha);

%%Step 2
P_2 = rotate_point_around_axis(P_1, Axis_Z_1_vector, Beta);
Axis_X_2_vector = rotate_axis_around_axis(Axis_X_1_vector,Axis_Z_1_vector,Beta);

%%Step 3

P_3 = rotate_point_around_axis(P_2, Axis_X_2_vector, Gamma);
end

