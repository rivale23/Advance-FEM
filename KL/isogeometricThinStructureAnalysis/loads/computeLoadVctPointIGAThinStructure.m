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
function Fl = computeLoadVctPointIGAThinStructure...
    (FlOld,xib,etab,p,q,Xi,Eta,CP,isNURBS,FAmp,direction,t,int,outMsg)
%% Function documentation
%
% Returns the consistent load vector for the application of a point load at
% the parametric location (xib,etab) to the Kirchhoff-Love shell.
%
%     Input :
%     FlOld : Existing force vector
%  xib,etab : Load position in the B-Spline parameter space
%       p,q : polynomial degree in u,v-direction
%    Xi,Eta : Knot vectors in xi,eta-direction
%        CP : Set of Control point coordinates and weights
%   isNURBS : Flag determining whether the basis is B-Spline or a NURBS
%      FAmp : The magnitude of the point load
% direction : The direction of FAmp, namely 1 = x, 2 = y, 3 = z
%         t : The time instance
%       int : On the integral integration (dummy variable for this 
%             function)
%    outMsg : Whether or not to output message
%             'outputEnabled' : enables output information   
%
%    Output :
%        Fl : The updated force vector
%
% Function layout :
%
% 0. Read input
%
% 1. Evaluate the NURBS basis functions at the point of the load application
%
% 2. Compute the load component at the application point
%
% 3. Re-assemble the load vector with respect to the global numbering
%
% 4. Update the existing load vector
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    if isvector(FlOld)
        fprintf('Update of the load vector corresponding to point\n');
    else
        fprintf('Computation of the load vector corresponding to point\n');
    end
    fprintf('load for the isogeometric Kirchhoff-Love shell problem\n');
    fprintf('has been initiated\n\n');
    fprintf('Xi parametric coordinate = %.2d\n',xib);
    fprintf('Eta parametric coordinate = %.2d\n',etab);
    fprintf('_______________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Find the knot span indices of the point where the load  is applied
xiSpan = findKnotSpan(xib,Xi,nxi);
etaSpan = findKnotSpan(etab,Eta,neta);

% Number of DOFs
noDOFs = 2*nxi*neta;

% Initialize output array
Fl = zeros(noDOFs,1);

%% 1. Evaluate the NURBS basis functions at the point of the load application
R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xib,Xi,etaSpan,q,etab,Eta,CP,isNURBS,0);

%% 2. Compute the load component at the application point

% Initialize the load array
FL = zeros(nxi,neta,3);

% Initialize counter
k = 1;

% Loop over all the contributions at the current knot span
for c = 0:q 
    for b = 0:p
        % Compute iteratively the load entry
        FL(xiSpan-p+b,etaSpan-q+c,direction) = R(k)*FAmp;
        
        % Update counter
        k = k + 1;
    end
end

%% 3. Re-assemble the load vector with respect to the global numbering
counter = 1;
for j = 1:length(FL(1,:,1))
    for i = 1:length(FL(:,1,1))
        % Assemble the x-coordinates of the load vector
        Fl(counter,1) = FL(i,j,1);
        
        % Assemble the y-coordinates of the load vector
        Fl(counter + 1,1) = FL(i,j,2);
        
        % Assemble the z-coordinates of the load vector
        Fl(counter + 2,1) = FL(i,j,3);
        
        % Update counter
        counter = counter + 3;
    end
end

%% 4. Update the existing load vector
if isvector(FlOld)  
    Fl = Fl + FlOld;   
end

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    if isvector(FlOld)
        fprintf('Load vector update took %.2d seconds \n\n',computationalTime);
        fprintf('________________Load Vector Update Ended_______________\n');
        fprintf('#######################################################\n\n\n');
    else
        fprintf('Load vector computation took %.2d seconds \n\n',computationalTime);
        fprintf('____________Load Vector Computation Ended______________\n');
        fprintf('#######################################################\n\n\n');
    end
end

end