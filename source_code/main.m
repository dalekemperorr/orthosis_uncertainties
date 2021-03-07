clc
clear

%% SELECT MODE
% 0- Perform calculations for positions read from file
% 1- Perform calculations for all possible positions
MODE = 0;


%% Nominal Orthosis axis definition [X, Y, Z, i, j, k] 
Axis_X_vector_nom = [0, 0, 0, 1, 0, 0];
Axis_Y_vector_nom = [0, 0, 0, 0, 1, 0];
Axis_Z_vector_nom = [0, 0, 0, 0, 0, 1];

%% Starting conditions:
Filename = 'Output_fromfile1.csv';
R = 150;
L = 350;
S = [0, 0, 0]';
R_v = [R, 0, 0]';
L_v = [0, 0, -L]';
A_kon = S + R_v;
B_kon = A_kon + L_v;
A_ort = A_kon; 
B_ort = B_kon;

%Perform calculations basing on selection
if(MODE == 0)
    [ Delta_Alpha, Delta_Beta, Delta_Gamma, Delta_T, GammaValues_vector ] = calculate_fromfile( Axis_X_vector_nom, Axis_Y_vector_nom, Axis_Z_vector_nom, A_kon, B_kon, A_ort, B_ort, S, R, L, Filename  );
else    
    [ Delta_Alpha, Delta_Beta, Delta_Gamma, Delta_T, GammaValues_vector ] = calculate_allpositions( Axis_X_vector_nom, Axis_Y_vector_nom, Axis_Z_vector_nom, A_kon, B_kon, A_ort, B_ort, S, R, L, Filename  );
end

extract_output2( Delta_Alpha, Delta_Beta, Delta_Gamma, Delta_T, GammaValues_vector );

disp('Calculations done, output data saved in file.')



