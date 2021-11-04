function run_example()
% Extract a resistance matrix from losses (different operating conditions).
%
%   Create a dummy resistance matrix.
%   Define different operating points (different current excitations).
%   Compute the associated losses (with the dummy resistance matrix).
%   Extract the resistance matrix from the operating points and the associated losses.
%   Compare the extracted resistance matrix with the dummy resistance matrix.
%
%   Thomas Guillod.
%   2021 - BSD License.

close('all')

%% Create a dummy resistance matrix for the different tests.

% A 3x3 matrix contains 6 independent coefficients.
R_11 = 1.6;
R_22 = 1.8;
R_33 = 2.0;
R_12 = 0.7;
R_13 = 0.5;
R_23 = 0.8;

% Construct the matrix (symmetric).
R_ref_mat = [R_11 R_12 R_13 ; R_12 R_22 R_23 ; R_13 R_23 R_33];

%% Extraction of the resistance matrix with a determined system.

% The excitation matrix defines the imposed currents (for different operating points).
% The excitation matrix contains 6 linearly independent operating points.
% Therefore, the equation system is determined.
I_operating_mat = [...
    1 0 0;...
    0 1 0;...
    0 0 1;...
    1 1 0;...
    0 1 1;...
    1 0 1;...
    ];

% Get the losses, extract the matrix, and compare the results with the reference matrix.
get_test('Determined system (6 equations, 6 variables, well-conditioned)', R_ref_mat, I_operating_mat);

%% Extraction of the resistance matrix with a quasi-singular system.

% The excitation matrix defines the imposed currents (for different operating points).
% The excitation matrix contains 6 linearly independent operating points.
% However, the different operating points are quasi-colinear.
% therefore, the equation system is determined but badly conditionned.
I_operating_mat = ones(10, 3)+1e-3.*rand(10,3);

% Get the losses, extract the matrix, and compare the results with the reference matrix.
get_test('Determined system (6 equations, 6 variables, ill-conditioned)', R_ref_mat, I_operating_mat);

%% Extraction of the resistance matrix with an over-determined system.

% The excitation matrix defines the imposed currents (for different operating points).
% The excitation matrix contains 10 linearly operating points.
% Therefore, the equation system is over-determined.
I_operating_mat = [...
    1 0 0;...
    0 1 0;...
    0 0 1;...
    1 1 0;...
    0 1 1;...
    1 0 1;...
    +1 -1 0;...
    0 +1 -1;...
    +1 0 -1;...
    1 1 1;...
    ];

% Get the losses, extract the matrix, and compare the results with the reference matrix.
get_test('Overdetermined system (10 equations, 6 variables)', R_ref_mat, I_operating_mat);

end

function get_test(tag, R_ref_mat, I_operating_mat)
% Test the resistance matrix extraction algorithm.
%
%    Parameters:
%        tag (string): description of the test 
%        R_ref_mat (matrix): reference resistance matrix 
%        I_operating_mat (matrix): considered operating conditions 

% Get the losses corresponding to the selection operating points.
P_operating_vec = 0.5.*diag(I_operating_mat*R_ref_mat*I_operating_mat');

% Extract the resistance matrix from the excitation matrix and the loss vector.
[R_operating_mat, res_vec, rcond_eqn] = get_res_ind_matrix(I_operating_mat, P_operating_vec);

% Scale the residuum into a relative residuum.
rel_res_vec = res_vec./P_operating_vec;

% Compute the relative error between the reference and extracted matrices.
rel_err_mat = (R_operating_mat-R_ref_mat)./R_ref_mat;

% Display the results.
fprintf('%s\n', tag)
get_disp(R_ref_mat, R_operating_mat, rcond_eqn, rel_res_vec, rel_err_mat);

end

function get_disp(R_ref_mat, R_operating_mat, rcond_eqn, rel_res_vec, rel_err_mat)
% Display the results of a resistance matrix extraction.
%
%    Parameters:
%        R_ref_mat (matrix): reference resistance matrix 
%        R_operating_mat (matrix): extracted resistance matrix 
%        rcond_eqn (scalar): equation system reciprocal condition 
%        rel_res_vec (vector): relative equation system residuum 
%        rel_err_mat (matrix): relative error between the reference and extracted matrices  

% Get the worst case errors.
rel_res_max = max(abs(rel_res_vec(:)));
rel_err_max = max(abs(rel_err_mat(:)));

% Compare the reference and extracted matrices.
fprintf('    Resistance matrix\n')
fprintf('        Reference resistance matrix = %s\n', mat2str(R_ref_mat, 5))
fprintf('        Extracted resistance matrix = %s\n', mat2str(R_operating_mat, 5))
fprintf('        Maximum relative error (resistance matrix coefficients) = %.3e\n', rel_err_max)

% Information of the equation system and the solution.
fprintf('    Equation system\n')
fprintf('        Equation system reciprocal condition = %.3e\n', rcond_eqn)
fprintf('        Maximum relative residuum (loss values) = %.3e\n', rel_res_max)

end
