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
function [x,hasConverged] = solve_LinearSystemMatlabBackslashSolver(A,b,x)
%
% Returns the solution to a system of linear equations using the standard
% matlab backslash operator
%
%   Input :
%       A : The left hand side matrix
%       b : The right hand side vector
%       x : The initial guess corresponding to an iterative solver (dummy 
%           for this function cause the matlab backslash operator 
%           corresponds to a direct solver)
%
%  Output :
%       x : The solution vector
%
% Function layout :
%
% 0. Read inout
%
% 1. Solve the linear system with the matlab backslash operator
%
% 2. Disable warning message
%
%% Function main body

%% 0. Read inout
hasConverged = true;
isWarningMessageEnabled = 'off';
% isWarningMessageEnabled = 'on';

%% 1. Solve the linear system with the matlab backslash operator
x = A\b;

%% 2. Disable warning message
% [msg, id] = lastwarn;
% if ischar(id)
%     if ~strcmp(id,'')
%         warning(isWarningMessageEnabled,id);
%     end
% else
%     warning(isWarningMessageEnabled,id);
% end

end