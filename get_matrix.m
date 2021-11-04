function [Q_mat, rcond_eqn, rel_res_vec] = get_matrix(I_mat, U_vec)

% get the number of equations
n_eqn = size(I_mat, 1);

% get the size of the matrix
n_mat = size(I_mat, 2);

% get the indices of the solution vector
[n_var, var_idx] = get_var_idx(n_mat);

% check size
assert(size(U_vec, 2)==1, 'parametrization is underdetermined')
assert(size(U_vec, 1)==n_eqn, 'parametrization is underdetermined')
assert(n_eqn>=n_var, 'parametrization is underdetermined')

% get the equations
for i=1:n_eqn
   eqn(i,:) = get_eqn(I_mat(i,:), var_idx);
end

% % check condition
if n_eqn==n_var
    rcond_eqn = rcond(eqn);
    assert(rcond_eqn>=eps, 'matrix is quasi-singular');
else
    rcond_eqn = min(svd(eqn))./max(svd(eqn));
    assert(rank(eqn)==n_var, 'matrix is quasi-singular');
end

% solve the equations
Q_vec = eqn\U_vec;

% assign the element to the matrix
Q_mat = get_mat(Q_vec, var_idx);

% get residum
U_sol_vec = 0.5.*diag(I_mat*Q_mat*I_mat');
rel_res_vec = (U_sol_vec-U_vec)./U_vec;

% get the min eigenvalue
eig_vec = eig(Q_mat);
assert(all(isreal(eig_vec)), 'unphysical complex eigenvalues')
assert(all(eig_vec>0), 'unphysical negatice eigenvalues')

end

function [n_var, var_idx] = get_var_idx(n_mat)

idx = 1;
for i=1:n_mat
    for j=1:i
        var_idx(idx, 1) = i;
        var_idx(idx, 2) = j;
        idx = idx+1;
    end
end
n_var = idx-1;

end

function eqn = get_eqn(I_vec, var_idx)

for i=1:size(var_idx, 1)
    idx_1 = var_idx(i,1);
    idx_2 = var_idx(i,2);
    
    if idx_1==idx_2
        eqn(:,i) = 0.5.*I_vec(idx_1).*I_vec(idx_2);
    else
        eqn(:,i) = 1.0.*I_vec(idx_1).*I_vec(idx_2);
    end
end

end

function Q_mat = get_mat(Q_vec, var_idx)

for i=1:size(var_idx, 1)
    idx_1 = var_idx(i,1);
    idx_2 = var_idx(i,2);
    
    Q_mat(idx_1, idx_2) = Q_vec(i);
    Q_mat(idx_2, idx_1) = Q_vec(i);
end

end