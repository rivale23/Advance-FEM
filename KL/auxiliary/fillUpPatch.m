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
function BSplinePatch = fillUpPatch(analysis,p,Xi,q,Eta,CP,isNURBS,...
    parameters,homDOFs,inhomDOFs,valuesInhomDOFs,weakDBC,cables,NBC,...
    masterDOFs,slaveDOFs,cb,xicoup,etacoup,int)
%% Function documentation
%
% Takes all nessecary arguments that a shell needs and assigns them to
% the structure.
%
%           Input :
%        analysis : .type : Analysis type
%             p,q : The polynomial degrees of the surface
%          Xi,Eta : The knot vectors in xi,eta-directions
%              CP : Set of Control Points and weights
%         isNURBS : Flag on whether the given patch is a NURBS or a 
%                   B-Spline patch
%      parameters : Technical ang geometrical information for the patch
%         homDOFs : Sequential numbering of the DOFs where homogeneous
%                   Dirichlet boundary conditions are applied
%       inhomDOFs : Sequential numbering of the inhomogeneous Dirichlet
%                   boundary conditions
% valuesInhomDOFs : Prescribed values on the inhomogeneous Dirichlet
%                  boundary conditions (in the same sequence as in the 
%                  array inhomDBC)
%         weakDBC : On the weak application of boundary conditions:
%                             .noCnd : Number of weak Dirichlet boundary
%                                      consitions
%                       .xiExtension : The extension of Dirichlet boundary
%                                      in the xi-direction
%                      .etaExtension : The extension of Dirichlet boundary
%                                      in the eta-direction
%                   .computeConstMtx : Function handle to the computation 
%                                      of the constant matrix accounting 
%                                      for the application of the weak 
%                                      boundary conditions
%              .computeTangMtxResVct : Function handle to the computation
%                                      of the tangent matrix and residual
%                                      vector accounting for the 
%                                      application of the weak boundary 
%                                      conditions
%          cables : Array related to cables that are embedded onto the
%                   B-Spline patch containing the following information :
%                           .noCables : Number of embedded cables
%                        .xiExtension : The xi-extension of the embedded 
%                                       cables in the B-Spline patch
%                       .etaExtension : The eta-extension of the embedded 
%                                       cables in the B-Spline patch
%                    .parametersCable : Array containing the technical and
%                                       geometrical parameters that a cable
%                                       must have :
%                                               .E : The Young's modulus of
%                                                    the cable
%                                        .radiusCS : The cross sectional
%                                                    radius of the cable
%                                             .rho : The density of the
%                                                    cable
%                                       .prestress : The prestress of the
%                                                    cable
%                                             .int : On the integration of
%                                                    the cables
%             NBC : On the Neumann boundary conditions:
%                             .noCnd : Number of Neumann boundary 
%                                      conditions
%                   .xiLoadExtension : Cell array {.noCnd} containing the 
%                                      load extensions in the xi-direction
%                  .etaLoadExtension : Cell array {.noCnd} containing the 
%                                      load extensions in the eta-direction
%                     .loadAmplitude : Array (1,.noCnd) containing the load
%                                      amplitudes
%                     .loadDirection : Array (1,.noCnd) containing the load
%                                      directions
%                    .computeLoadVct : Cell array {.noCnd} containing the 
%                                      function name for the computation of
%                                      the load vector
%                    .isConservative : Array (1,.noCnd) of flags indicating
%                                      whether the load is conservative or
%                                      not
%      masterDOFs : Array containing the global numbering of the DOFs which
%                   are master in master-slave relation. Used when for
%                   example DOF i is forced to be equal to DOF j
%       slaveDOFs : Array containing the global numbering of the DOFs which
%                   are forced to be equal to the masterDOFs in the given
%                   sequence, that is, slaveDOFs[i] = slaveDOFs[i]
%              cb : DOFs on coupling interface
%  xicoup,etacoup : Arrays containing all coupled regions
%             int : On the numerical integration
%
%          Output :
%    BSplinePatch : Structure containing all above information together
%                   with the following:
%                              .noCPs : Number of Control Points
%                           .noElmnts : Number of elements
%                         .xiKnotSpan : Array containing the xi- knot span
%                                       indices for each element
%                        .etaKnotSpan : Array containing the xi- knot span
%                                       indices for each element
%                          .xiIndexCP : The index of the Control points in
%                                       the xi-direction that affect each
%                                       element
%                         .etaIndexCP : The index of the Control points in
%                                       the eta-direction that affect each
%                                       element
%                            .noGPsEl : Number of Gauss Points per element
%                                 .GW : Array containing all Gauus Point 
%                                       weights, namely, xiGW*etaGW for a 
%                                       surface patch
%                    .elementAreaOnGP : Array containing the product of the
%                                       determinant of the transformation 
%                                       from the parent to the parameter 
%                                       space with the determinant of the 
%                                       transformation from the parameter 
%                                       space to the physical space and the
%                                       Gauss point weight at each Gauss
%                                       point
%                            .xi,.eta : Array containing the images of the 
%                                       Gauss Points on the parameter space
%                                  .R : Array containing the basis
%                                       functions at the element level for
%                                       each Gauss point and at each 
%                                       element
%                      .dRdXi,.dRdEta : Arrays containing the first 
%                                       derivatives of the basis functions 
%                                       at the element level for each Gauss 
%                                       point and at each element
%                          .GXi,.GEta : Arrays containing the values of the
%                                       base vectors along each parametric
%                                       direction for each Gauss point
%               .prestressVoigtVector : Array containing the prestress 
%                                       values for each element and Gauss 
%                                       point
%                            .G3Tilde : Array containing the surface normal
%                                       (not normalized vector) for each
%                                       Gauss Point
%                          .minElArea : The minimum element area for the
%                                       isogeometric discretization for the
%                                       given surface patch
%                          .maxElArea : The maximum element area for the
%                                       isogeometric discretization for the
%                                       given surface patch
%                                .EFT : Array containing the Element
%                                       Freedom Tables for each element of
%                                       the  isogeometric discretization 
%                                       for the given surface patch
%                               .p,.q : Polynomial orders of the surface
%                                       patch
%                            .Xi,.Eta : The knot vectors of the surface
%                                       patch
%                                 .CP : The set of Control Point
%                                       coordinates and weights for the
%                                       surface patch
%                                .CPd : The set of the displaced Control 
%                                       Point coordinates and weights for 
%                                       the surface patch
%                            .isNURBS : Flag on whether the basis is a
%                                       B-Spline or a NURBS
%                         .parameters : The technical parameters of the
%                                       surface patch
%                            .homDOFs : The global numbering of the
%                                       homogeneous Dirichlet boundary
%                                       conditions
%                          .inhomDOFs : The global numbering of the
%                                       inhomogeneous Dirichlet boundary
%                                       conditions
%                    .valuesInhomDOFs : The values of the prescribed DOFs
%                                       at the DOFs where inhomogeneous
%                                       Dirichlet boundary conditions are
%                                       applied
%                            .weakDBC : On the application of weak boundary
%                                       conditions
%                             .cables : On the cables that are embedded in
%                                       the B-Spline
%                                .NBC : Array containing the information on
%                                       the Neumann boundary conditions
%                         .masterDOFs : The master DOFs
%                          .slaveDOFs : The slave DOFs
%                                 .cb : Array containing the global
%                                       numbering of the DOFs which lie on
%                                       a coupling boundary
%                                .int : Array containing information on the
%                                       integration
%                    .xicoup,.etacoup : Coupling information along the xi
%                                       and the eta parameter lines
%                       .DOFNumbering : Array containing the numbering of 
%                                       the DOFs
%                        .boundingBox : The bounding box containing the
%                                       B-Spline patch
%
% Function layout :
%
% 0. Read input
%
% 1. Create a DOF numbering for the B-Spline patch
%
% 2. Compute the bounding box of the B-Spline geometry
%
% 3. Get the Gauss point coordinates and weights according the selected quadrature rule
% ->
%    3i. Save the span index of the current element and the indices of the Control Points affecting the current knot span
%
%   3ii. Create the EFT of the current element
%
%  3iii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
%
%   3iv. Initialize element area and counter for the Gauss Points in the element level
%
%    3v. Loop over all Gauss points
%    ->
%        3v.1. Compute and save the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%        3v.2. Compute and save the NURBS basis functions and their derivatives according to the chosen analysis
%
%        3v.3. Compute and save the covariant base vectors of the reference configuration
%
%        3v.4. Compute and save the surface normal of the reference configuration (third covariant base vector not normalized)
%
%        3v.5. Compute and save the Gauss weight
%
%        3v.6 Compute the element area on the Gauss Point and add the contribution
%
%        3v.7. Update the GP counter        
%    <-
%   3vi. Find the minimum element area in the isogemetric mesh
%
%  3vii. Update element counter
% <-
%
%% Function main body

%% 0. Read input

% Assign NURBS parameters and boundary conditions to the patch array
BSplinePatch.p = p;
BSplinePatch.q = q;
BSplinePatch.Xi = Xi;
BSplinePatch.Eta = Eta;
BSplinePatch.CP = CP;
BSplinePatch.CPd = BSplinePatch.CP;
BSplinePatch.isNURBS = isNURBS;
BSplinePatch.parameters = parameters;
BSplinePatch.homDOFs = homDOFs;
BSplinePatch.inhomDOFs = inhomDOFs;
BSplinePatch.valuesInhomDOFs = valuesInhomDOFs;
BSplinePatch.weakDBC = weakDBC;
BSplinePatch.cables = cables;
BSplinePatch.NBC = NBC;
BSplinePatch.masterDOFs = masterDOFs;
BSplinePatch.slaveDOFs = slaveDOFs;
BSplinePatch.cb = cb;
BSplinePatch.int = int;
BSplinePatch.xicoup = xicoup;
BSplinePatch.etacoup = etacoup;

% Issue warnings for analyses which have not yet been optimized
if strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis') || ...
        strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
    warning('Analysis type %s has not yet been optimized',analysis.type);
end

% Number of knots in xi- and eta- directions
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in xi- and eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check the input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Number of Control Points to the patch array
BSplinePatch.noCPs = nxi*neta;

% Number of Degrees of Freedom (DOFs) per element
if strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis') || ...
        strcmp(analysis.type,'isogeometricMembraneAnalysis')
    noDOFsNode = 3;
elseif strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
    noDOFsNode = 2;
end
noCPsEl = (p+1)*(q+1);
noDOFsEl = noDOFsNode*noCPsEl;

% Issue the Gauss Point coordinates and weights
if isfield(int,'type')
    if strcmp(int.type,'default')
        if strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis') || ...
            strcmp(analysis.type,'isogeometricMembraneAnalysis') || ...
            strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
            noXiGP = p + 1;
            noEtaGP = q + 1;
        end
    elseif strcmp(int.type,'user')
        noXiGP = int.xiNGP;
        noEtaGP = int.etaNGP;
    end
else
    return;
end
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(noXiGP);
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(noEtaGP);

% Number of Gauss Points per element
BSplinePatch.noGPsEl = noXiGP*noEtaGP;

% Get the number of the elements in the isogeometric discretization
noElmntsXi = length(unique(Xi)) - 1;
noElmntsEta = length(unique(Eta)) - 1;
BSplinePatch.noElmnts = noElmntsXi*noElmntsEta;

% Initialize the span indices for each element
BSplinePatch.xiKnotSpan = zeros(BSplinePatch.noElmnts,1);
BSplinePatch.etaKnotSpan = zeros(BSplinePatch.noElmnts,1);

% Total number of Gauss Points
noGPs = BSplinePatch.noGPsEl*BSplinePatch.noElmnts;

% Initialize the element area size on the Gauss point for each Gauss point
BSplinePatch.elementAreaOnGP = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl);

% Initialize arrays of Gauss point coordinates, basis functions, their
% derivatives, base vectors and their derivatives
BSplinePatch.xi = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl);
BSplinePatch.eta = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl);
BSplinePatch.R = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl,noCPsEl);
if strcmp(analysis.type,'isogeometricMembraneAnalysis')
    BSplinePatch.dRdXi = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl,noCPsEl);
    BSplinePatch.dRdEta = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl,noCPsEl);
    BSplinePatch.GXi = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl,3);
    BSplinePatch.GEta = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl,3);
    BSplinePatch.G3Tilde = zeros(3,noGPs);
    BSplinePatch.prestressPage = zeros(BSplinePatch.noElmnts,BSplinePatch.noGPsEl,3);
end

% Initialize minimum element area in the IGA mesh
tolerance = 1e-4;
if abs(CP(1,1,1)-CP(nxi,1,1)) >= tolerance
    BSplinePatch.minElArea = abs(CP(1,1,1)-CP(nxi,1,1));
elseif abs(CP(1,1,1)-CP(1,neta,1)) >= tolerance
    BSplinePatch.minElArea = abs(CP(1,1,1)-CP(1,neta,1));
elseif abs(CP(1,end,1)-CP(nxi,end,1))
    BSplinePatch.minElArea = abs(CP(1,end,1)-CP(nxi,end,1));
else
    BSplinePatch.minElArea = abs(CP(end,1,1)-CP(end,neta,1));
end
BSplinePatch.maxElArea = 0;

% Initialize counters
counterElmnts = 1;
counterGPs = 1;

% Initialize the array of the Element Freedom Tables (EFTs)
BSplinePatch.EFT = zeros(noDOFsEl,BSplinePatch.noElmnts);

%% 1. Create a DOF numbering for the B-Spline patch
if strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis') || ...
        strcmp(analysis.type,'isogeometricMembraneAnalysis')
    BSplinePatch.DOFNumbering = zeros(nxi,neta,3);
    k = 1;
    for cpj = 1:neta
        for cpi = 1:nxi
            BSplinePatch.DOFNumbering(cpi,cpj,1) = k;
            BSplinePatch.DOFNumbering(cpi,cpj,2) = k + 1;
            BSplinePatch.DOFNumbering(cpi,cpj,3) = k + 2;
            k = k + 3;
        end
    end
elseif strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
    BSplinePatch.DOFNumbering = zeros(nxi,neta,2);
    k = 1;
    for cpj = 1:neta
        for cpi = 1:nxi
            BSplinePatch.DOFNumbering(cpi,cpj,1) = k;
            BSplinePatch.DOFNumbering(cpi,cpj,2) = k + 1;
            k = k + 2;
        end
    end
end

%% 2. Compute the bounding box of the B-Spline geometry
BSplinePatch.boundingBox = [BSplinePatch.CP(1,1,1)
                            BSplinePatch.CP(1,1,1)
                            BSplinePatch.CP(1,1,2)
                            BSplinePatch.CP(1,1,2)
                            BSplinePatch.CP(1,1,3)
                            BSplinePatch.CP(1,1,3)];
for counterXi = 1:length(BSplinePatch.CP(:,1,1))
    for counterEta = 1:length(BSplinePatch.CP(1,:,1))
        x = BSplinePatch.CP(counterXi,counterEta,1);
        y = BSplinePatch.CP(counterXi,counterEta,2);
        z = BSplinePatch.CP(counterXi,counterEta,3);
        if x < BSplinePatch.boundingBox(1,1)
            BSplinePatch.boundingBox(1,1) = x;
        elseif x > BSplinePatch.boundingBox(2,1)
            BSplinePatch.boundingBox(2,1) = x;
        end
        if y < BSplinePatch.boundingBox(3,1)
            BSplinePatch.boundingBox(3,1) = y;
        elseif y > BSplinePatch.boundingBox(4,1)
            BSplinePatch.boundingBox(4,1) = y;
        end
        if z < BSplinePatch.boundingBox(5,1)
            BSplinePatch.boundingBox(5,1) = z;
        elseif z > BSplinePatch.boundingBox(6,1)
            BSplinePatch.boundingBox(6,1) = z;
        end
    end
end

%% 3. Get the Gauss point coordinates and weights according the selected quadrature rule
for iEtaSpan = q+1:meta-q-1
    for iXiSpan = p+1:mxi-p-1
        % check if element is greater than zero
        if Xi(iXiSpan+1) ~= Xi(iXiSpan) && Eta(iEtaSpan+1) ~= Eta(iEtaSpan)
            %% 3i. Save the span index of the current element and the indices of the Control Points affecting the current knot span
            BSplinePatch.xiKnotSpan(counterElmnts,1) = iXiSpan;
            BSplinePatch.etaKnotSpan(counterElmnts,1) = iEtaSpan;
            BSplinePatch.xiIndexCP(counterElmnts,:) = ....
                BSplinePatch.xiKnotSpan(counterElmnts,1) - p + (0:p);
            BSplinePatch.etaIndexCP(counterElmnts,:) = ....
                BSplinePatch.etaKnotSpan(counterElmnts,1) - q + (0:q);
            
            %% 3ii. Create the EFT of the current element
            if strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis') || ...
                strcmp(analysis.type,'isogeometricMembraneAnalysis')
                k = 1;
                for cpj = iEtaSpan-q:iEtaSpan
                    for cpi = iXiSpan-p:iXiSpan
                        BSplinePatch.EFT(k,counterElmnts) = BSplinePatch.DOFNumbering(cpi,cpj,1);
                        BSplinePatch.EFT(k+1,counterElmnts) = BSplinePatch.DOFNumbering(cpi,cpj,2);
                        BSplinePatch.EFT(k+2,counterElmnts) = BSplinePatch.DOFNumbering(cpi,cpj,3);
                        k = k + 3;
                    end
                end
            elseif strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
                k = 1;
                for cpj = iEtaSpan-q:iEtaSpan
                    for cpi = iXiSpan-p:iXiSpan
                        BSplinePatch.EFT(k,counterElmnts) = BSplinePatch.DOFNumbering(cpi,cpj,1);
                        BSplinePatch.EFT(k+1,counterElmnts) = BSplinePatch.DOFNumbering(cpi,cpj,2);
                        k = k + 2;
                    end
                end
            end
            
            %% 3iii. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
            %
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(iXiSpan+1)-Xi(iXiSpan))*(Eta(iEtaSpan+1)-Eta(iEtaSpan))/4;
            
            %% 3iv. Initialize element area and counter for the Gauss Points in the element level
            elementArea = 0;
            counterGPEl = 1;
            
            %% 3v. Loop over all Gauss points
            for cEta = 1:noEtaGP
                for cXi = 1:noXiGP
                    %% 3v.1. Compute and save the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                    xi = ( Xi(iXiSpan+1)+Xi(iXiSpan) + xiGP(cXi)*(Xi(iXiSpan+1)-Xi(iXiSpan)) )/2;
                    eta = ( Eta(iEtaSpan+1)+Eta(iEtaSpan) + etaGP(cEta)*(Eta(iEtaSpan+1)-Eta(iEtaSpan)) )/2;
                    BSplinePatch.xi(counterElmnts,counterGPEl) = xi;
                    BSplinePatch.eta(counterElmnts,counterGPEl) = eta;
                    
                    %% 3v.2. Compute and save the NURBS basis functions and their derivatives according to the chosen analysis
                    if strcmp(analysis.type,'isogeometricMembraneAnalysis') || ...
                            strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
                        noDrvBasis = 1;
                    elseif strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis')
                        noDrvBasis = 2;
                    end
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (iXiSpan,p,xi,Xi,iEtaSpan,q,eta,Eta,CP,isNURBS,noDrvBasis);
                    BSplinePatch.R(counterElmnts,counterGPEl,:) = dR(:,1);
                    if strcmp(analysis.type,'isogeometricMembraneAnalysis')
                        BSplinePatch.dRdXi(counterElmnts,counterGPEl,:) = dR(:,2);
                        BSplinePatch.dRdEta(counterElmnts,counterGPEl,:) = dR(:,3);
                    end
                    
                    %% 3v.3. Compute and save the covariant base vectors of the reference configuration
                    if strcmp(analysis.type,'isogeometricMembraneAnalysis') || ...
                            strcmp(analysis.type,'isogeometricPlateInMembraneActionAnalysis')
                        nDrvBaseVct = 0;
                    elseif strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis')
                        nDrvBaseVct = 1;
                    end
                    [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (iXiSpan,p,iEtaSpan,q,CP,nDrvBaseVct,dR);
                    if strcmp(analysis.type,'isogeometricMembraneAnalysis')
                        BSplinePatch.GXi(counterElmnts,counterGPEl,:) = dG1(:,1);
                        BSplinePatch.GEta(counterElmnts,counterGPEl,:) = dG2(:,1);
                    end
                    
                    %% 3v.4. Compute and save the prestress values
                    if strcmp(analysis.type,'isogeometricMembraneAnalysis')
                        isPrestressOverDefinedSystem = false;
                        if (isfield(parameters.prestress,'computeParametricCoordinates') && ~isfield(parameters.prestress,'computeBaseVectors')) || ...
                                (~isfield(parameters.prestress,'computeParametricCoordinates') && isfield(parameters.prestress,'computeBaseVectors'))
                            error('Both or none of function handles parameters.prestress.computeParametricCoordinates and parameters.prestress.computeBaseVectors have to be defined');
                        elseif isfield(parameters.prestress,'computeParametricCoordinates') && isfield(parameters.prestress,'computeBaseVectors')
                            isPrestressOverDefinedSystem = true;
                        end
                        if isPrestressOverDefinedSystem
                            GabCovariant = [dG1(:,1) dG2(:,1)]'*[dG1(:,1) dG2(:,1)];
                            GContravariant = GabCovariant\[dG1(:,1) dG2(:,1)]';
                            GContravariant = GContravariant';
                            eLC = computeLocalCartesianBasis4BSplineSurface([dG1(:,1) dG2(:,1)],GContravariant);
                            X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                                (iXiSpan,p,xi,Xi,iEtaSpan,q,eta,Eta,CP,dR(:,1));
                            theta = parameters.prestress.computeParametricCoordinates(X);
                            prestressBaseVct = parameters.prestress.computeBaseVectors(theta(1,1),theta(2,1));
                            T2LC = computeT2LocalCartesianBasis(prestressBaseVct,eLC);
                            BSplinePatch.prestressVoigtVector(counterElmnts,counterGPEl,:) = ...
                                T2LC*parameters.prestress.voigtVector;
                        else
                            BSplinePatch.prestressVoigtVector(counterElmnts,counterGPEl,:) = ...
                                parameters.prestress.voigtVector;
                        end
                    end
                 
                    %% 3v.4. Compute and save the surface normal of the reference configuration (third covariant base vector not normalized)
                    BSplinePatch.G3Tilde(:,counterGPs) = cross(dG1(:,1),dG2(:,1));
                    dA = norm(BSplinePatch.G3Tilde(:,counterGPs));
                    
                    %% 3v.5. Compute and save the Gauss weight
                    GW = xiGW(cXi)*etaGW(cEta);
                    
                    %% 3v.6 Compute the element area on the Gauss Point and add the contribution
                    BSplinePatch.elementAreaOnGP(counterElmnts,counterGPEl) = dA*detJxiu*GW;
                    elementArea = elementArea + BSplinePatch.elementAreaOnGP(counterElmnts,counterGPEl);
                    
                    %% 3v.7. Update the GP counter
                    counterGPEl = counterGPEl + 1;
                    counterGPs = counterGPs + 1;
                end
            end
            %% 3vi. Find the minimum element area in the isogemetric mesh
            if elementArea < BSplinePatch.minElArea
                BSplinePatch.minElArea = elementArea;
            end
            if elementArea > BSplinePatch.maxElArea
                BSplinePatch.maxElArea = elementArea;
            end
            
            %% 3vii. Update element counter
            counterElmnts = counterElmnts + 1;
        end
    end
end

end