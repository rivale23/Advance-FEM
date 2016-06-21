%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,stiffMtx,stiffMtx_disturbed,dindex] = solve_IGALinearSystem_shortcut...
    (BSplinePatch,Position,u,...
    freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
    solve_LinearSystem)
%% Function documentation
%
% Returns the solution to a linear system which correspond to the
% isogeometric discretization of the underlying field.
%
% Function layout :
%
% 1. Compute the linear matrices of the system
%
% 2. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
%
% 3. Solve the linear equation system
%
% 4. Re-assemble to the complete vector of unknowns
%
%% Function main body

%% 1. compute stiffness matrix and disturbed stiffness matrix
[stiffMtx,stiffMtx_disturbed,RHS,dindex,a] = computeLinearMtrcsSensitivity(BSplinePatch,Position);

%% 2. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
if norm(valuesInhomDOFs(inhomDOFs)) ~= 0
    RHS = RHS - stiffMtx(:,inhomDOFs)*valuesInhomDOFs(inhomDOFs);
end

%% 3. Solve the linear equation system
[uRed,hasLinearSystemConverged] = solve_LinearSystem...
    (stiffMtx(freeDOFs,freeDOFs),RHS(freeDOFs),u(freeDOFs));
if ~hasLinearSystemConverged
    error('Linear equation solver has not converged');
end

%% 4. Re-assemble to the complete vector of unknowns
u(freeDOFs) = uRed;
u(homDOFs) = 0;
u(inhomDOFs) = valuesInhomDOFs(inhomDOFs);

end