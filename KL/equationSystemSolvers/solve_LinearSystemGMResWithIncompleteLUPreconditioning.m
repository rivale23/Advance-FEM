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
function [x,hasConverged] = solve_LinearSystemGMResWithIncompleteLUPreconditioning...
    (A,b,x)
%% Function documentation
%
% Returns the solution to a system of linear equations Ax=b using the GMRES 
% (Generalized Minimal Residual) solver with an incomplete LU factorization
%
%	Input :
%       A : The left hand side matrix
%       b : The right hand side vector
%       x : The initial guess of the GMRes iterations
%
%  Output :
%       x : The solution vector
%
%   Function layout :
%
% 0. Read input
%
% 1. Perform the LU factorization
%
% 2. Loop over the GMRes iterations
% ->
%    2i. Solve the GMRes step
%
%   2ii. Check convergence criterion
% <-
% 
% 3. Check if the GMRes solver conveged
%
%% Function main body

%% 0. Read input

% set the parameters of the iterations
tolerance = 1e-12;
maxIt = 100;

% Initialize the iteration counter
counterIter = 1;

% Initialize convergence flag
hasConverged = true;

%% 1. Perform the LU factorization
[L,U] = ilu(sparse(A));

%% 2. Loop over the GMRes iterations
while counterIter < maxIt
    %% 2i. Solve the GMRes step
    [x,~,relResidual] = gmres(A,b,[],tolerance,[],L,U,x);
    counterIter = counterIter + 1;
    
    %% 2ii. Check convergence criterion
    if relResidual < tolerance
        break;
    end
end

%% 3. Check if the GMRes solver conveged
if counterIter == maxIt;
    warning( 'GMRes solver did not converge in the given number of iterations!' );
    hasConverged = false;
end

end