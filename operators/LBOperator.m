function S = LBOperator(mesh, basis)
% Returns the functional representation of the LB operator
% let [v_1,...,v_n] = basis
% f = sum_i ( a_i v_i )
% Lf = sum_i ( a_i L v_i ) = sum_i sum_j (a_j <L v_j, v_i> v_i) = 
%    = sum_i b_i v_i
% where b_i = sum_j (a_j <L v_j, v_i>)
% So, b = S a , S_ij = <L v_j, v_i>, 
% or if not orthogonal: L v_j = sum_i (S_ij v_i)

S = descriptorsToCoeffs(LB(mesh)*basis, basis);
S(abs(S)<1e-8)=0;
