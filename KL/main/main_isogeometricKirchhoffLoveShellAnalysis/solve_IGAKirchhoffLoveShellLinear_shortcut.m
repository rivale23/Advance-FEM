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
function [dHat,stiffMtx,stiffMtx_disturbed,dindex] = solve_IGAKirchhoffLoveShellLinear_shortcut...
    (BSplinePatch,Position,solve_LinearSystem)
%% Function documentation
% 
% Returns the displacement field and the complete force vector of a
% linear Kirchhoff-Love shell. Source Reference : 
%
% J. Kiendl "Isogeometric Analysis and Shape Optimal Design for Shell 
% Structures" Ph.D. Thesis, Technische Universtät München (2011)
%
% Function layout:
%
% 0. Read input
%
% 1. Solve the linear system
%
%% Function main body

%% 0. Read input

t = 0;

% Re-assign the arrays
CP = BSplinePatch.CP;

% Number of Control Points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Create an element freedom table for the patch in the array
BSplinePatch.DOFNumbering = zeros(nxi,neta,3);
k = 1;
for cpj = 1:neta
    for cpi = 1:nxi
        BSplinePatch.DOFNumbering(cpi,cpj,1) = k;
        BSplinePatch.DOFNumbering(cpi,cpj,2) = k + 1;
        BSplinePatch.DOFNumbering(cpi,cpj,3) = k + 2;

        % Update counter
        k = k + 3;
    end
end

% Create the element freedom table for the BSplinePatch into the array of
% the patches
BSplinePatch.EFTPatches = 1:3*BSplinePatch.noCPs;

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
dHat = zeros(noDOFs,1);

%% 1. Solve the linear system
[dHat,stiffMtx,stiffMtx_disturbed,dindex] = solve_IGALinearSystem_shortcut...
    (BSplinePatch,Position,dHat,...
    freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
    solve_LinearSystem,t);

end
