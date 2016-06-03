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
function Fl = computeLoadVctAreaIGAThinStructure...
    (FlOld,xib,etab,p,q,Xi,Eta,CP,isNURBS,FAmp,direction,t,int,outMsg)
%
% Returns the consistent nodal forces to an area load fload [N/m^2] for the 
% isogeometric Kirchhoff-Love shell problem. The direction of f can be in 
% x, y, z, parallel or perpendicular to the edge.
%
%     Input :
%     FlOld : Existing force vector
%  xib,etab : load extension (e.g. xib = [0 1], etab = 1)
%       p,q : the polynomial degrees of the surface
%    Xi,Eta : the knot vectors of the surface
%        CP : The Control Point coordinates and weights
%   isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
%      FAmp : The amplitude of the constant line load or handle to load 
%             function [N/m] for varying loads
% direction : Direction of the applied force:
%               1-x, 2-y, 3-z, self-weight
%               4-x, 5-y, 6-z, snow load
%               7-pressure load (perpendicular to the shell)
%         t : The time instance
%       int : On the numerical integration
%    outMsg : Whether or not to output message on refinement progress
%             'outputEnabled' : enables output information
%
%    Output :
%        Fl : Updated force vector
%
% Function layout :
%
% 0. Read input
%
% 1. On the numerical integration
%
% 2. Find the parametrization of the area integral
%
% 3. Loop over the elements
% ->
%    3i. Compute the determinant of the Jacobian to the transformation from the parameter to the integration domain
%
%   3ii. Loop over all Gauss points
%   ->
%        3ii.1. Compute the coordinates, the map and the Gauss Point location and weight for the fixed parametric coordinate
%
%        3ii.2. Compute the IGA basis functions and their first derivatives
%
%        3ii.3. Compute the base vectors of the configuration
%
%        3ii.4. the physical location of the load application if it is applied a varying load
%
%        3ii.5. Compute the product R(xi,eta)*dX in the physical space at the Gauss Point
%
%        3ii.6. Compute the element load vector at the quadrature point
%
%        3ii.7. Assemble the element load vector to the global load vector
%   <-
% <-
%
% 4. Re-assemble the load vector with respect to the global numbering
%
% 5. Update the existing load vector
%
% 6. Appendix
% 
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    if isvector(FlOld)
        fprintf('Update of the load vector corresponding to area\n');
    else
        fprintf('Computation of the load vector corresponding to area\n');
    end
    fprintf('load for the isogeometric Kirchhoff-Love shell problem\n');
    fprintf('has been initiated\n\n');
    if isnumeric(FAmp)
        fprintf('Constant area load is assumed with amplitude = %.2d\n',FAmp);
    else
        fprintf('Varying area load is assumed\n');
    end
    if isscalar(xib)
        fprintf('Load extension in xi parametric direction = [%.2d,%.2d]\n',xib,xib);
    else
        fprintf('Load extension in xi parametric direction = [%.2d,%.2d]\n',xib(1),xib(2));
    end
    if isscalar(etab)
        fprintf('Load extension in eta parametric direction = [%.2d,%.2d]\n',etab,etab);
    else
        fprintf('Load extension in eta parametric direction = [%.2d,%.2d]\n',etab(1),etab(2));
    end
    fprintf('_______________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Length of the knot vectors
mxi = length(Xi);
meta = length(Eta);

% Number of control points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% if the fload parameter is numeric then assign directly its value
if isnumeric(FAmp) == 1
    FAmplitude = FAmp;  
end

% Number of DOFs
nDOFs = 2*nxi*neta;

% Initialize auxiliary array
Rdx = zeros(p+1,q+1);

% Initialize global force vector
F = zeros(nxi,neta,3);

% Initialize output array
Fl = zeros(nDOFs,1);

%% 1. On the numerical integration

% Initialize structure
nGP = zeros(2,1);

% Issue the Gauss points for the selected integration scheme:
if strcmp(int.type,'default')
   % Default scheme is the full gaussian quadrature element-wise (FGI)
   nGP(1) = ceil((p+1)/2);
   nGP(2) = ceil((q+1)/2);
elseif strcmp(int.type,'manual')
    % Manual choice of the gauss points
    nGP(1) = int.xiNGPForLoad;
    nGP(2) = int.xetaNGPForLoad;
end

% Issue the Gauss points for the numerical integration
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(nGP(1));
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(nGP(2));

%% 2. Find the parametrization of the area integral

% For the xi-direction :
% ______________________

% Start and end knot spans in xi-direction for the load application
i1 = findKnotSpan(xib(1),Xi,nxi);
i2 = findKnotSpan(xib(2),Xi,nxi);

% If the end span of the integration domain is not the end of the knot
% span decrease the knot span index by 1 (I have no clue why Kiendl
% implemented that in this way)
if xib(2)~=Xi(mxi)  
    i2 = i2 - 1;  
end

% For the eta-direction :
% _______________________

% Start and end knot spans in eta-direction for the load application
j1 = findKnotSpan(etab(1),Eta,neta);
j2 = findKnotSpan(etab(2),Eta,neta);

% If the end span of the integration domain is not the end of the knot
% span decrease the knot span index by 1 (I have no clue why Kiendl
% implemented that in this way)
if etab(2)~=Eta(meta)
    j2 = j2 - 1;  
end
                
%% 3. Loop over the elements
for j = j1:j2
    for i = i1:i2
        % check if element is greater than zero
        if Xi(i+1)~=Xi(i) && Eta(j+1)~=Eta(j)
            %% 3i. Compute the determinant of the Jacobian to the transformation from the parameter to the integration domain
            detJParam2Integr = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;
            
            %% 3ii. Loop over all Gauss points
            for cEta = 1:nGP(2)
                for cXi = 1:nGP(1)
                 %% 3ii.1. Compute the coordinates, the map and the Gauss Point location and weight for the fixed parametric coordinate
                    
                    % Transform the xi-coordinate to the intergation domain
                    xi = ( Xi(i+1)+Xi(i) + xiGP(cXi)*(Xi(i+1)-Xi(i)) )/2;
                    
                    % Transform the v-coordinate to the intergation domain
                    eta = ( Eta(j+1)+Eta(j) + etaGP(cEta)*(Eta(j+1)-Eta(j)) )/2;
                    
                    % The Gauss weight of the 2-dimensional quadrature
                    GW = xiGW(cXi)*etaGW(cEta);
                    
                  %% 3ii.2. Compute the IGA basis functions and their first derivatives
                    nDrv = 1;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,nDrv);
                    
                    %% 3ii.3. Compute the base vectors of the configuration
                    
                    % Compute the in-plane base vectors
                    nDrv = 0;
                    [G1,G2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i,p,j,q,CP,nDrv,dR);
                    
                    % Compute the not normalized surface normal
                    N = cross(G1,G2);
                    
                    % Compute the differential element area
                    A = norm(N);
                    
                    % Compute the surface normal
                    n = N/A;
                    
                    %% 3ii.4. the physical location of the load application if it is applied a varying load
                    if ~isnumeric(FAmp)
                        % Compute the Cartesian coordinates on the B-Spline
                        % surface
                        S = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                            (i,p,xi,Xi,j,q,eta,Eta,CP,dR(:,1));
                        x = S(1,1);
                        y = S(1,1);
                        z = S(1,1);
                        
                        % Compute the load amplitude on the Gauss Point is the load is not constant
                        if ischar(FAmp) && ~isa(FAmp,'function_handle')
                            FAmp = str2func(FAmp);
                        end
                        FAmplitude = FAmp(x,y,z,t);
                    else
                        FAmplitude = FAmp;
                    end
                    
                    %% 3ii.5. Compute the product R(xi,eta)*dX in the physical space at the Gauss Point
                    
                    % Initialize counter
                    counter = 1;
                    
                    % Compute its values recursively
                    for c = 0:q
                        for b = 0:p
                            Rdx(b+1,c+1) = dR(counter,1)*A;

                            % update counter
                            counter = counter + 1;
                        end
                    end
                    
                    %% 3ii.6. Compute the element load vector at the quadrature point
                    Fel = FAmplitude*detJParam2Integr*GW*Rdx;
                    
                    %% 3ii.7. Assemble the element load vector to the global load vector
                    if direction==1 || direction==2 || direction==3
                        F(i-p:i,j-q:j,direction) = Fel(:,:) + F(i-p:i,j-q:j,direction);
                    elseif direction==4 || direction==5 || direction==6
                        F(i-p:i,j-q:j,direction-3) = Fel(:,:)*n(direction-3) + F(i-p:i,j-q:j,direction-3);
                    elseif direction==7
                        for idir = 1:3
                            F(i-p:i,j-q:j,idir) = Fel(:,:)*n(idir) + F(i-p:i,j-q:j,idir);
                        end
                    end
                end
            end
        end 
    end
end

%% 4. Re-assemble the load vector with respect to the global numbering
counter = 1;
for j = 1:length(F(1,:,1))
    for i = 1:length(F(:,1,1))
        % Assemble the x-coordinates of the load vector
        Fl(counter,1) = F(i,j,1);
        
        % Assemble the y-coordinates of the load vector
        Fl(counter + 1,1) = F(i,j,2);
        
        % Assemble the z-coordinates of the load vector
        Fl(counter + 2,1) = F(i,j,3);
        
        % Update counter
        counter = counter + 3;
    end
end

%% 5. Update the existing load vector
if isvector(FlOld)  
    Fl = Fl + FlOld;   
end

%% 6. Appendix
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