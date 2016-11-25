function [s r] = collectOperators(mesh1, mesh2, options)
% Collects operator commutativity constraints from mesh1 and mesh2.  
% Returns operators as nxn matrices

s = {};
r = {};

no = length(options.operators);
for i=1:no
    oper = options.operators{i};
    
    s{i} = getOp(mesh1, oper);
    r{i} = getOp(mesh2, oper);
end

function s = getOp(mesh, name)        
if strcmp(name, 'LB')
    s = LB(mesh);

elseif strcmp(name, 'confLB')
    s = mesh.cotLaplace;

elseif strcmp(name, 'SymGT')
    load([mesh.name '_sym']);
    s = sparse(1:mesh.nv,corr,ones(mesh.nv,1),mesh.nv,mesh.nv);
    
else fprintf(1,'Unrecognized operator: %s', descriptor);
end
