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
function [Xir,Etar,CPr] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,nxi,neta,outMsg)
%% Function documentation
%
% Computes the CP's and U after equidistant knot insertion into the given
% curve
%
%    Input :
% nxi,neta : Number of knots to be inserted equidistantly in xi-,eta 
%            direction
%      p,q : The polynomial degrees of the surface
%   Xi,Eta : The knot vectors of the NURBS surface in xi,eta-direction
%       CP : The Control Points of the NURBS surface
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
%   Output :
%      Xir : The knot vector in xi-direction after the refinement
%     Etar : The knot vector in eta-direction after the refinement
%      CPr : The set of Control points after the refinement
%
% Function layout :
%
% 0. Read input
%
% 1. Create the vectors with the knots to be added in both parametric directions
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________________________\n');
    fprintf('################################################################\n');
    fprintf('Uniform Knot insertion for a B-Spline curve has been initiated \n\n');
    fprintf('Number of knots before knot insertion in xi-direction nxi = %d\n',length(Xi));
    fprintf('Number of knots after knot insertion in xi-direction nxi = %d\n',length(Xi)+nxi);
    fprintf('Number of knots before knot insertion in eta-direction neta = %d\n',length(Eta));
    fprintf('Number of knots after knot insertion in eta-direction neta = %d\n',length(Eta)+neta);
    fprintf('________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Initialize the vectors containing the additional knots 
Rxi = zeros(nxi-1);
Reta = zeros(neta-1);

%% 1. Create the vectors with the knots to be added in both parametric directions

% Loop and add knots into the vector in u-direction
for i = 1:nxi-1  
    Rxi(i) = i/nxi*(Xi(length(Xi))-Xi(1));  
end

% Loop and add knots into the vector in v-direction
for i = 1:neta-1  
    Reta(i) = i/neta*(Eta(length(Eta))-Eta(1));  
end

% Apply knot insertion onto the surface
[Xir,Etar,CPr] = knotRefineBSplineSurface(p,Xi,q,Eta,CP,Rxi,Reta,'');

%% 2. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Knot insertion took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Knot Insertion Ended______________________\n');
    fprintf('################################################################\n\n\n');
end

end