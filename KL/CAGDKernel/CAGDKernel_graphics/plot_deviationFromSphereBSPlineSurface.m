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
function [relErrorInL2,index] = plot_deviationFromSphereBSPlineSurface...
    (BSplinePatches,dHat,propSphere,int,graph,outMsg)
%% Function documentation
%
% Returns the relative error in the L2-norm of the difference of the given 
% B-Spline surface from an exact sphere as well as the index of the current 
% graph.
%
%           Input :
%  BSplinePatches : Array of B-Spline patches each of which contains :
%                       .p,.q : The polynomial orders of the B-Spline patch
%                    .Xi,.Eta : The knot vectors of the B-Spline patch
%                         .CP : The set of Control point coordinates and
%                               weights of the B-Spline patch
%                 .EFTPatches : The freedom tables for each B-Spline patch
%                               in the multipatch geometry
%            dHat : The displacement field corresponding to the multipatch
%                   geometry, if chosen as a string the variable CPd of the
%                   B-Spline patch is chosen
%      propSphere : Analytical data for the analytical description of the
%                   sphere :
%                       .center : The center of the exact sphere
%                       .radius : The radius of the exact sphere
%             int : Quadrature scheme for the computation of the L2-norm :
%                       .type : 'default' or 'user'
%                     .noGPXi : Number of Gauss points for the integration
%                               in the xi-direction if int.type = 'user'
%                    .noGPEta : Number of Gauss points for the integration
%                               in the eta-direction if int.type = 'user'
%           graph : On the graphics
%          outMsg : Enables outputting information onto the command window
%                   when chosen as 'outputEnabled'
%
%          Output :
%    relErrorInL2 : The relative error in the L2-norm from the deviation of
%                   the shape from a sphere
%           index : The index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all patches in the B-Spline patch geometry
% ->
%    1i. Get the patch parameters
%
%   1ii. Compute the displaced Control Points of the patch
%
%  1iii. Get the Gauss points and weights of the quadrature for the computation of the L2-norm
%
%   1iv. Initialize the array of the deviation from a sphere for the plot
%
%    1v. Loop over all the elements
%    ->
%        1v.1. Create an element freedom table
%
%        1v.2. Loop over all Gauss points
%
%        1v.3. Loop over all Gauss points
%        ->
%              1v.3i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%             1v.3ii. Compute the NURBS basis functions and their first derivatives at the Gauss Point
%
%            1v.3iii. Compute the covariant base vectors of the reference configuration at the Gauss point
%
%             1v.3iv. Compute the surface normal of the reference configuration at the Gauss point (third covariant base vector not normalized)
%
%              1v.3v. Compute the legth of G3Tilde at the Gauss point (= area dA of the undeformed configuration)
%
%             1v.3vi. Compute the point on the deformed B-Spline surface
%
%            1v.3vii. Compute the direction vector of the computed point on the B-Spline surface
%
%           1v.3viii. Find the corresponding spherical coordinates with respect to the computed direction vector
%
%             1v.3ix. Compute the coordinates of the point on the sphere
%
%              1v.3x. Compute the element area on the Gauss Point
%
%             1v.3xi. Compute the L2-norm of the difference from a sphere and the norm of the exact sphere
%        <-
%    <-
%
%   1vi. Initialize the eta parametric coordinate and the counter in eta-direction
%
%  1vii. Loop over all the sampling points in -eta direction
%  ->
%        1vii.1. Find the span in the eta-direction
%
%        1vii.2. Initialize the xi parametric coordinate and the counter in xi-direction
%
%        1vii.3. Loop over all the sampling points in -xi direction
%        ->
%                1vii.3i. Find the span in xi-direction
%
%               1vii.3ii. Compute the IGA basis functions and their derivatives
%
%              1vii.3iii. Compute the point on the deformed B-Spline surface
%
%               1vii.3iv. Compute the direction vector of the computed point on the B-Spline surface
%
%                1vii.3v. Find the corresponding spherical coordinates with respect to the computed direction vector
%
%              1vii.3vii. Compute the deviation from the exact sphere at the samping point
%
%             1vii.3viii. Update the parametric coordinate and the counter in xi-direction
%        <-
%
%        1vii.4. pdate the parametric coordinate in eta-direction and the counter    
%  <-
%
% 1viii. Plot the deformed surface together with the difference from an exact sphere
% <-
%
% 2. Assign graphic properties
%
% 3. Compute the relative error in the L2-norm for the deviation from the exact sphere
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_________________________________________________________________________\n');
    fprintf('#########################################################################\n');
    fprintf('Plotting of the distribution of the deviation from an exact sphere over a\n');
    fprintf('multipatch B-Spline surface and computation of the corresponding L2-norm \n');
    fprintf('has been initiated \n');
    fprintf('_________________________________________________________________________\n\n');
end
%% 0. Read input

% Initialize variables
errorL2 = 0;
exactL2 = 0;

% Number of points for the visualization of the deviation from the sphere
xiGrid = 49;
etaGrid = 49;

% Assign a tolerance value
tol = 10e-10;

% Number of patches
noPatches = length(BSplinePatches);

% Initialize handle to the figure
figure(graph.index)

%% 1. Loop over all patches in the B-Spline patch geometry
for iPatches = 1:noPatches
    %% 1i. Get the patch parameters
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    isNURBS = BSplinePatches{iPatches}.isNURBS;
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    mxi = length(Xi);
    mxiUnique = length(unique(Xi));
    meta = length(Eta);
    metaUnique = length(unique(Eta));
    noDOFsEl = 3*(p + 1)*(q + 1);
    DOFNumbering = BSplinePatches{iPatches}.DOFNumbering;
    EFTPatches = BSplinePatches{iPatches}.EFTPatches;
    dxi = (Xi(mxi)-Xi(1))/(etaGrid-1);
    deta = (Eta(meta)-Eta(1))/(xiGrid-1);
    
    %% 1ii. Compute the displaced Control Points of the patch
    if ~ischar(dHat)
        CPd = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
            (CP,dHat(EFTPatches));
    else
        CPd = BSplinePatches{iPatches}.CPd;
    end
    
    %% 1iii. Get the Gauss points and weights of the quadrature for the computation of the L2-norm
    if strcmp(int.type,'default')
        noGPXi = ceil((p + 1)/2);
        noGPEta = ceil((q + 1)/2);
        isQuadratureUser = false;
    elseif strcmp(int.type,'user')
        noGPXi = int.noGPXi;
        noGPEta = int.noGPEta;
        [GPXi,GWXi] = getGaussPointsAndWeightsOverUnitDomain(noGPXi);
        [GPEta,GWEta] = getGaussPointsAndWeightsOverUnitDomain(noGPEta);
        isQuadratureUser = true;
    else
        error('Define the type of quadrature rule');
    end
    if ~isQuadratureUser
        [GPXi,GWXi] = getGaussPointsAndWeightsOverUnitDomain(noGPXi);
        [GPEta,GWEta] = getGaussPointsAndWeightsOverUnitDomain(noGPEta);
    end
    
	%% 1iv. Initialize the array of the deviation from a sphere for the plot
    clear deviationFromSphere X;
    deviationFromSphere = zeros(noGPXi*mxiUnique,noGPEta*metaUnique);
    P = zeros(noGPXi*mxiUnique,noGPEta*metaUnique,3);
    
    %% 1v. Loop over all the elements
    for j = q+1:meta-q-1
        if Eta(j+1) ~= Eta(j)
            for i = p+1:mxi-p-1
                if Xi(i+1) ~= Xi(i)
                    %% 1v.1. Create an element freedom table

                    % Initialize element freedome table
                    EFT = zeros(1,noDOFsEl);

                    % initialize counter
                    k = 1;

                    % relation global-local dof
                    for cpj = j-q:j
                        for cpi = i-p:i
                            EFT(k) = DOFNumbering(cpi,cpj,1);
                            EFT(k+1) = DOFNumbering(cpi,cpj,2);
                            EFT(k+2) = DOFNumbering(cpi,cpj,3);

                            % Update counter
                            k = k + 3;
                        end
                    end

                    %% 1v.2. Loop over all Gauss points
                    %
                    %         | xi_i+1 - xi_i                    |
                    %         | -------------            0       |
                    %         |        2                         |
                    %  xi,u = |                                  |
                    %         |                  eta_j+1 - eta_j |
                    %         |        0         --------------- |
                    %         |                          2       |
                    detJxiu = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;

                    %% 1v.3. Loop over all Gauss points
                    for iGPEta = 1:noGPEta
                        for iGPXi = 1:noGPXi
                            %% 1v.3i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                            xi = ( Xi(i+1)+Xi(i) + GPXi(iGPXi)*(Xi(i+1)-Xi(i)) )/2;
                            eta = ( Eta(j+1)+Eta(j) + GPEta(iGPEta)*(Eta(j+1)-Eta(j)) )/2;

                            %% 1v.3ii. Compute the NURBS basis functions and their first derivatives at the Gauss Point
                            dR = computeIGABasisFunctionsAndDerivativesForSurface...
                                (i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,1);

                            %% 1v.3iii. Compute the covariant base vectors of the reference configuration at the Gauss point
                            [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                                (i,p,j,q,CP,0,dR);

                            %% 1v.3iv. Compute the surface normal of the reference configuration at the Gauss point (third covariant base vector not normalized)
                            A3Tilde = cross(A1(:,1),A2(:,1));

                            %% 1v.3v. Compute the legth of G3Tilde at the Gauss point (= area dA of the undeformed configuration)
                            dA = norm(A3Tilde);

                            %% 1v.3vi. Compute the point on the deformed B-Spline surface
                            X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                                (i,p,xi,Xi,j,q,eta,Eta,CPd,dR(:,1));

                            %% 1v.3vii. Compute the direction vector of the computed point on the B-Spline surface
                            directionVct = X - propSphere.center;

                            %% 1v.3viii. Find the corresponding spherical coordinates with respect to the computed direction vector
                            theta = atan(directionVct(2,1)/directionVct(1,1));
                            phi = acos(directionVct(3,1)/norm(directionVct));

                            %% 1v.3ix. Compute the coordinates of the point on the sphere
                            XSphere = [propSphere.radius*cos(theta)*sin(phi)
                                       propSphere.radius*sin(theta)*sin(phi)
                                       propSphere.radius*cos(phi)];

                            %% 1v.3x. Compute the element area on the Gauss Point
                            elementAreaOnGP = dA*detJxiu*GWXi(iGPXi)*GWEta(iGPEta);

                            %% 1v.3xi. Compute the L2-norm of the difference from a sphere and the norm of the exact sphere
                            errorL2 = errorL2 + norm(abs(X) - abs(XSphere))^2*elementAreaOnGP;
                            exactL2 = exactL2 + norm(XSphere)^2*elementAreaOnGP;
                        end
                    end
                end
            end
        end
    end
    
    %% 1vi. Initialize the eta parametric coordinate and the counter in eta-direction
    eta = Eta(1);
    etaCounter = 1;  
    
    %% 1vii. Loop over all the sampling points in -eta direction
    while eta <= Eta(meta)+tol
        %% 1vii.1. Find the span in the eta-direction
        etaSpan = findKnotSpan(eta,Eta,neta);
        
        %% 1vii.2. Initialize the xi parametric coordinate and the counter in xi-direction
        xi = Xi(1);
        xiCounter = 1;
        
        %% 1vii.3. Loop over all the sampling points in -xi direction
         while xi <= Xi(mxi)+tol
            %% 1vii.3i. Find the span in xi-direction
            xiSpan = findKnotSpan(xi,Xi,nxi);
            
            %% 1vii.3ii. Compute the IGA basis functions and their derivatives
            R = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
            
            %% 1vii.3iii. Compute the point on the deformed B-Spline surface
            P(xiCounter,etaCounter,1:3) = ...
                computeCartesianCoordinatesOfAPointOnBSplineSurface...
                (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CPd,R);
            
            %% 1vii.3iv. Compute the direction vector of the computed point on the B-Spline surface
            directionVct = squeeze(P(xiCounter,etaCounter,:)) - propSphere.center;
            
            %% 1vii.3v. Find the corresponding spherical coordinates with respect to the computed direction vector
            theta = atan(directionVct(2,1)/directionVct(1,1));
            phi = acos(directionVct(3,1)/norm(directionVct));
             
            %% 1vii.3vi. Compute the coordinates of the point on the sphere
            XSphere = [propSphere.radius*cos(theta)*sin(phi)
                       propSphere.radius*sin(theta)*sin(phi)
                       propSphere.radius*cos(phi)];
                   
            %% 1vii.3vii. Compute the deviation from the exact sphere at the samping point
            unitVct = (XSphere - propSphere.center)/norm(XSphere - propSphere.center);
            deviationFromSphere(xiCounter,etaCounter) = ...
                (abs(squeeze(P(xiCounter,etaCounter,:))) - abs(XSphere))'*unitVct;
            
            %% 1vii.3viii. Update the parametric coordinate and the counter in xi-direction
            xi = xi + dxi;
            xiCounter = xiCounter + 1;
         end
         
        %% 1vii.4. pdate the parametric coordinate in eta-direction and the counter
        eta = eta + deta;
        etaCounter = etaCounter + 1;
    end
    
    %% 1viii. Plot the deformed surface together with the difference from an exact sphere
    surf(P(:,:,1),P(:,:,2),P(:,:,3),deviationFromSphere(:,:));
    hold on;
    plot_knotsForBSplineSurfaceOnCartesianSpace...
        (p,q,Xi,Eta,CPd,isNURBS,false,xiGrid,etaGrid);
end
hold off;

%% 2. Assign graphic properties
shading interp;
title('Distribution of the deviation from an exact sphere');
colormap('jet');
colorbar;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
index = graph.index + 1;

%% 3. Compute the relative error in the L2-norm for the deviation from the exact sphere
relErrorInL2 = sqrt(errorL2)/sqrt(exactL2);

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Plotting and computation took %.2d seconds \n\n',computationalTime);
    fprintf('__________________Plotting Current Configuration Ended___________________\n');
    fprintf('#########################################################################\n\n\n');
end

end