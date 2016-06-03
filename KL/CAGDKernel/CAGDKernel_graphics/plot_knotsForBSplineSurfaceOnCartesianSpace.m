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
function plot_knotsForBSplineSurfaceOnCartesianSpace...
    (p,q,Xi,Eta,CP,isNURBS,isDeformed,xiGrid,etaGrid)
%% Function documentation
%
% Draws the element edges for the NURBS surface, i.e the knots on the
% geometry
%
%      Input :
%        p,q : Polynomial degrees
%     Xi,Eta : Knot vectors in u,v-direction
%         CP : Control Point coordinates and weights
%    isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
% isDeformed : Flag on whether the plot is called for the reference or the
%              current configuration
%     xiGrid : Points to use in xi-direction
%    etaGrid : Points to use in eta-direction
%
%     Output : 
%              Graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the edges in xi-direction
%
% 2. Plot the edges in eta-direction
%
% 3. Plot all the element edges
%
%% Function main body

%% 0. Read input

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in u,v-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Assign a tolerance value
eps = 10e-10;

% Initialize counter
l = 1;

% Compute step size in xi-direction
dxi = (Xi(mxi)-Xi(1))/xiGrid;

% Compute step size in eta-direction
deta = (Eta(meta)-Eta(1))/etaGrid; 

% Initialize plotting array
P = zeros(xiGrid,etaGrid,3);

%% 1. Plot the edges in xi-direction
for j2 = q+1:meta-q
    % Get the starting coordinate in eta-direction
    eta = Eta(j2);
    
    % Find the span in eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    % Get the starting coordinate in xi-direction
    xi = Xi(1);
    
    % Initialize counter for the edges in xi-direction
    k = 1;
    
    % Loop over all the coordinates in xi-direction
    while xi <= Xi(mxi) + eps
        % Find the span in xi-direction
        xiSpan = findKnotSpan(xi,Xi,nxi);
        
        % Compute the IGA basis functions in xi-direction
        nDrv = 0;
        R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        % Compute the Cartesian image of the parametric point (xi,eta)
        P(k,l,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
        
        % Update counter for the edges in xi-direction
        k = k + 1;
        
        % Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    % Update the counter for the edges in eta-direction
    l = l + 1;
end

%% 2. Plot the edges in eta-direction
for i2 = p+1:mxi-p
    % Get the starting coordinate in xi-direction
    xi = Xi(i2);
    
    % Find the span in xi-direction
    xiSpan = findKnotSpan(xi,Xi,nxi);
    
    % Get the starting coordinate in eta-direction
    eta = Eta(1);
    
    % Initialize counter for the edges in eta-direction
    k = 1;
    
    % Loop over all the coordinates in eta-direction
    while eta <= Eta(meta)+eps
        % Find the span in eta-direction
        etaSpan = findKnotSpan(eta,Eta,neta);
        
        % Compute the IGA basis functions in xi-direction
        nDrv = 0;
        R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        % Compute the Cartesian image of the parametric point (xi,eta)
        P(k,l,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
        
        % Update counter for the edges in eta-direction
        k = k + 1;
        
        % Update the parametric coordinate in eta-direction
        eta = eta + deta;
    end
    % Update the counter for the edges in xi-direction
    l = l + 1;
end

%% 3. Plot all the element edges
if isDeformed == 0
    plot3(P(:,:,1),P(:,:,2),P(:,:,3),'Color','black','LineWidth',.01);
    axis equal;
    grid on;
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    zlabel('z','FontSize',18);
else
    plot3(P(:,:,1),P(:,:,2),P(:,:,3),'Color','black','LineStyle','-.');
    axis equal;
    grid on;
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    zlabel('z','FontSize',18);
end

end