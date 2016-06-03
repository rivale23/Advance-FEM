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
function [X,Y,Z] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,xiGrid,etaGrid)
%% Function documentation
%
% Returns three arrays, containing the Cartesian coordinates of the points 
% on the NURBS surface in a grid of xiGrid Times etaGridv lines
%
%   Input :
%     p,q : Polynomial degrees
%  Xi,Eta : Knot vectors in xi,eta-direction
%      CP : Control point coordinates and weights
% isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
%  xiGrid : Number of lines in xi-direction
% etaGrid : Number of lines in eta-direction
%
%  Output :
%       X : Array containing the x-coordinates of the points on the surface
%       Y : Array containing the y-coordinates of the points on the surface
%       Z : Array containing the z-coordinates of the points on the surface
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the parametric coordinate locations
% ->
%     1i. Compute start coordinate on the xi parameter lines and initialize counter
%
%    1ii. Compute the knot span index in eta-direction
%
%   1iii. Loop over all the coordinates in xi-direction
%   ->
%         1iii.1.  Compute the knot span index in xi-direction
%
%         1iii.2. Compute the IGA basis functions at (xi,eta)
%
%         1iii.3. Compute the Cartesian coordinates of (xi,eta)
%
%         1iii.4. Update counter for the lines in xi-direction
%
%         1iii.5. Update the xi-parametric coordinate
%   <-
%    1iv. Update counter for the lines in eta-direction
%
%     1v. Update the eta-parametric coordinate
% <-
% 2. Write the coordinates into the individual arrays
%
%% Function main body

%% 0. Read input

% Number of knots in xi,eta-direction
mxi = length(Xi);
meta = length(Eta);

% Number of control points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check the compatibility of the NURBS parameters
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Assign a tolerance value
eps = 10e-10;

% Initialize counter for the lines in eta-direction
lEta = 1;  

% Incremental step for the lines in eta-direction
deta = (Eta(meta)-Eta(1))/etaGrid;    

% Incremental step for the lines in xi-direction
dxi = (Xi(mxi)-Xi(1))/xiGrid;  

% Start coordinate on the eta parameter lines
eta = Eta(1);

% Initialize output array
S = zeros(xiGrid,etaGrid,3);

%% 1. Loop over all the parametric coordinate locations
while eta <= Eta(meta)+eps
    %% 1i. Compute start coordinate on the xi parameter lines and initialize counter
    xi = Xi(1);
    
    % Initialize counter for the lines in u-direction
    lxi = 1;
    
    %% 1ii. Compute the knot span index in eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    %% 1iii. Loop over all the coordinates in xi-direction
    while xi <= Xi(mxi)+eps
        %% 1iii.1.  Compute the knot span index in xi-direction
    	xiSpan = findKnotSpan(xi,Xi,nxi);
    	
        %% 1iii.2. Compute the IGA basis functions at (xi,eta)
        nDrv = 0;
        R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        %% 1iii.3. Compute the Cartesian coordinates of (xi,eta)
        S(lxi,lEta,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
    	
        %% 1iii.4. Update counter for the lines in xi-direction
    	lxi = lxi + 1;
        
        %% 1iii.5. Update the xi-parametric coordinate
        xi = xi + dxi;
    end
  %% 1iv. Update counter for the lines in eta-direction
  lEta = lEta + 1;
  
  %% 1v. Update the eta-parametric coordinate
  eta = eta + deta;
end

%% 2. Write the coordinates into the individual arrays
X = S(:,:,1);
Y = S(:,:,2);
Z = S(:,:,3);

end