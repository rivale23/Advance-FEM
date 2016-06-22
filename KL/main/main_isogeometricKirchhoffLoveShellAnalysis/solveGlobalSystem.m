function [ u ] = solveGlobalSystem( K_global, F_global, BSplinePatch, solve_LinearSystem)
%SOLVEGLOBALSYSTEM Summary of this function goes here
%   Detailed explanation goes here

%% 0. Read input

% Re-assign the arrays
CP = BSplinePatch.CP;

% Number of Control Points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Get number of DOFs
noDOFs = 3*nxi*neta;

% Find the numbering of the DOFs where homogeneous Dirichlet conditions are
% prescribed
homDOFs = BSplinePatch.homDOFs;

% Find the numbering of the free DOFs
freeDOFs = zeros(noDOFs,1);
for i=1:noDOFs
    freeDOFs(i,1) = i;
end
freeDOFs(ismember(freeDOFs,homDOFs)) = [];

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Initialize the displacement field
u = zeros(noDOFs,1);

%% 1. Solve the linear equation system
[uRed,hasLinearSystemConverged] = solve_LinearSystem...
    (K_global(freeDOFs,freeDOFs),F_global(freeDOFs),u(freeDOFs));
if ~hasLinearSystemConverged
    error('Linear equation solver has not converged');
end

%% 2. Re-assemble to the complete vector of unknowns
u(freeDOFs) = uRed;
u(homDOFs) = 0;
u(inhomDOFs) = valuesInhomDOFs(inhomDOFs);

end

