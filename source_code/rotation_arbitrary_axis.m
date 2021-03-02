%==========================================================================
%[name] rotate_point_around_axis
%[desc] Rotate input point around axis that is defined by two arbitrary points. 
%       The rotation is defined by input angle by
%[in]   P1 - veector [x y z] describing input point's position 
%[in]   Pa0 -  vector[x y z] describing first arbitrary point located on axis of rotation
%[in]   Pa1 -  vector[x y z] describing second arbitrary point located on axis of rotation
%[in]   fi - angle of rotation in radians
%[out]  P2 - vector[x y z ] describing rotated point's orientation
%==========================================================================
function P2 =rotation_arbitrary_axis(P1, Pa0, Pa1, fi)

P2 = zeros(length(P1), 1); 

A = Pa1(1) - Pa0(1);  
B = Pa1(2) - Pa0(2);  
C = Pa1(3) - Pa0(3); 

L = sqrt(A^2 + B^2 + C^2);
V = sqrt(B^2 + C^2);

D = eye(4,4);
D(1,4) = -Pa0(1);
D(2,4) = -Pa0(2);
D(3,4) = -Pa0(3);

Rx = eye(4,4);
Rx(2,2) = C/V;
Rx(2,3) = -B/V;
Rx(3,2) = B/V;
Rx(3,3) = C/V;

Ry = eye(4,4);
Ry(1,1) = V/L;
Ry(1,3) = -A/L;
Ry(3,1) = A/L;
Ry(3,3) = V/L;

Rz = eye(4,4);
Rz(1,1) = cos(fi);
Rz(1,2) = -sin(fi);
Rz(2,1) = sin(fi);
Rz(2,2) = cos(fi);

T = inv(D) * inv(Rx) * inv(Ry) * Rz * Ry * Rx * D;
P1_ext = ones(4,1);
P1_ext(1) = P1(1);
P1_ext(2) = P1(2);
P1_ext(3) = P1(3);


P2_ext = T * P1_ext;

P2(1) = P2_ext(1);
P2(2) = P2_ext(2);
P2(3) = P2_ext(3);
  
end
