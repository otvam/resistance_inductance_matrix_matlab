function [Q_mat, res_vec, rcond_eqn] = get_matrix(I_mat, U_vec)
% Extract a resistance/inductance matrix from operating points.
%
%    Different operating points are provided (current excitations).
%    The associated losses/energies are provided.
%    The corresponding resistance/inductance matrix is extracted.
%
%    The following the quadratic form is used
%        - U = 0.5*I'*Q_mat*I
%        - I is the current excitation vector
%        - U is the loss/energy value
%        - Q_mat is the resistance/inductance matrix
%
%    The following numbers describe the number of variables:
%        - n_mat: size of the resistance/inductance matrix
%        - n_var: number of independent resistance/inductance matrix coefficients
%        - n_op: number of provided operating points
%        - n_var = (n_mat*(n_mat+1))/2
%
%    The size of the different matrices and vectors are:
%        - I_mat (n_op x n_mat): current excitation matrix
%        - U_vec (n_op x 1): loss/energy vector
%        - Q_mat (n_mat x n_mat): resistance/inductance matrix
%        - res_vec (n_op x 1): residuum vector for the loss/energy
%
%    The resistance/inductance matrix is extracted with the following method:
%        - n_op<n_var: under-determined equation system => invalid problem
%        - n_op==n_var: determined equation system => linear equation system
%        - n_op>n_var: over-determined equation system => least-square fit
%
%    Errors are raised for the following cases:
%        - the equation system is under-determined 
%        - the equation system singular or quasi-singular
%        - the resistance/inductance matrix has negative eigenvalues
%
%    Parameters:
%        I_mat (matrix): current excitation matrix 
%        U_vec (vector): loss/energy vector
%
%    Returns:
%        Q_mat (matrix): resistance/inductance matrix
%        res_vec (vector): residuum vector for the loss/energy
%        rcond_eqn (scalar): equation system reciprocal condition
%
%    Thomas Guillod.
%    2021 - BSD License.

% get the number operating points
n_op = size(I_mat, 1);

% get the size of the resistance/inductance matrix
n_mat = size(I_mat, 2);

% get the number of independent coefficient and their indices
[n_var, var_idx] = get_var_idx(n_mat);

% check data size
assert(size(U_vec, 2)==1, 'parametrization is underdetermined')
assert(size(U_vec, 1)==n_op, 'parametrization is underdetermined')
assert(n_op>=n_var, 'parametrization is underdetermined')

% get the equations for the different operating conditions
for i=1:n_op
   eqn(i,:) = get_eqn(I_mat(i,:), var_idx);
end

% get the reciprocal condition (svd decomposition)
rcond_eqn = min(svd(eqn))./max(svd(eqn));

% check if the equation system is not quasi-singular
assert(rank(eqn)==n_var, 'matrix is quasi-singular');
assert(rcond_eqn>eps, 'matrix is quasi-singular');

% solve the equation system
Q_vec = eqn\U_vec;

% assign the solution vector to the resistance/inductance matrix
Q_mat = get_mat(Q_vec, var_idx);

% check the the eigenvalues are positive
assert(all(eig(Q_mat)>0), 'unphysical negative eigenvalues')

% compute the residum
U_sol_vec = 0.5.*diag(I_mat*Q_mat*I_mat');
res_vec = U_sol_vec-U_vec;

end

function [n_var, var_idx] = get_var_idx(n_mat)
% Get the number of independent coefficient and their indices.
%
%    The resistance/inductance matrix (n_mat x n_mat) is symmetric.
%    The resistance/inductance matrix contains n_var independent coefficient.
%
%    Parameters:
%        n_mat (scalar): size of the inductance matrix 
%
%    Returns:
%        n_var (scalar): number of independent coefficient in the resistance/inductance matrix
%        var_idx (matrix): matrix tracking the resistance/inductance matrix indices

% number of independent coefficient in the resistance/inductance matrix
n_var = (n_mat.*(n_mat+1))./2;

% creating a matrix mapping the following indices:
%     - the linear variable indices: 1:n_var
%     - the position in the resistance/inductance matrix: i,j
idx = 1;
for i=1:n_mat
    for j=1:i
        var_idx(idx, 1) = i;
        var_idx(idx, 2) = j;
        idx = idx+1;
    end
end

end

function eqn = get_eqn(I_vec, var_idx)
% Get the equation corresponding to a specific operating point.
%
%    Parameters:
%        I_vec (vector): current vector for the considered operating point
%        var_idx (matrix): matrix tracking the resistance/inductance matrix indices
%
%    Returns:
%        eqn (vector): vector corresponding to a line of the equation system

% for each independent coefficient in the resistance/inductance matrix
for i=1:size(var_idx, 1)
    % get the indices of this coefficient in the resistance/inductance matrix
    idx_1 = var_idx(i,1);
    idx_2 = var_idx(i,2);
    
    % assign the coefficient of the equation system
    if idx_1==idx_2
        eqn(:,i) = 0.5.*I_vec(idx_1).*I_vec(idx_2);
    else
        eqn(:,i) = 1.0.*I_vec(idx_1).*I_vec(idx_2);
    end
end

end

function Q_mat = get_mat(Q_vec, var_idx)
% Assign the solution vector to the resistance/inductance matrix.
%
%    Parameters:
%        Q_vec (vector): solution vector with the coefficients of the resistance/inductance matrix
%        var_idx (matrix): matrix tracking the resistance/inductance matrix indices
%
%    Returns:
%        Q_mat (matrix): resistance/inductance matrix

% for each independent coefficient in the resistance/inductance matrix
for i=1:size(var_idx, 1)
    % get the indices of this coefficient in the resistance/inductance matrix
    idx_1 = var_idx(i,1);
    idx_2 = var_idx(i,2);
    
    % assign the coefficient and the symmetric counterpart
    Q_mat(idx_1, idx_2) = Q_vec(i);
    Q_mat(idx_2, idx_1) = Q_vec(i);
end

end