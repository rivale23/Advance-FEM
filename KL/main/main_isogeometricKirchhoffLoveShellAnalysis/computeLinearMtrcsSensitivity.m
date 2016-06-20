
%takes 
function [K,KDist,F,EFTDist,minElArea] = computeLinearMtrcsSensitivity(BSplinePatch,disturbed_cp)

if isempty(disturbed_cp)
    disturbed_cp = [1,1];
end

%% 0. Read input

% Assign back the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
CPDist= BSplinePatch.CPd;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
int = BSplinePatch.int;
DOFNumbering = BSplinePatch.DOFNumbering;

% Neuman boundary conditions
NBC = BSplinePatch.NBC;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
% Number of control points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Initialize minimum element area in the IGA mesh
tolerance = 1e-4;
if abs(CP(1,1,1)-CP(nxi,1,1))>=tolerance
    minElArea = abs(CP(1,1,1)-CP(nxi,1,1));
else
    minElArea = CP(1,1,1)-CP(1,neta,1);
end

% Local number of DOFs
noDOFsEl = 3*(p+1)*(q+1);
% Global number of DOFs
noDOFs = 3*nxi*neta;

% Initialize global stiffness matrix
K  = zeros(noDOFs,noDOFs);
% Initialize global disturbed stiffness matrix
KDist  = zeros(noDOFs,noDOFs);

%% 1. Compute the material matrices

% Compute the membrane material matrix
Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
    [1              parameters.nue 0
	 parameters.nue 1              0
     0              0              (1-parameters.nue)/2];
                                                 
% Compute the bending material matrix
Db = parameters.E*parameters.t^3/(12*(1-parameters.nue^2))*...
    [1              parameters.nue 0
     parameters.nue 1              0 
     0              0              (1-parameters.nue)/2];

%% 2. Choose an integration rule

% Select the integration scheme
if strcmp(int.type,'default')
    xiNGP = p + 1;
    etaNGP = q + 1;
elseif strcmp(int.type,'manual')
    xiNGP = int.xiNGP;
    etaNGP = int.etaNGP;
end

% Issue the Gauss Point coordinates and weights
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(xiNGP);
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(etaNGP);

%% 3. loops over elements

% check which elements will be affected by distortion
% span of indices of basis functions overlapping with the changed basis function (i,j)
[cpi_span, cpj_span] = meshgrid( max(1,disturbed_cp(1)-p) : min(mxi-p-1,disturbed_cp(1)+p), max(1,disturbed_cp(2)-q) : min(meta-q-1,disturbed_cp(2)+q));
% initialize corresponding element freedom table                     
EFTDist = zeros(1,numel(cpi_span)*3);
 
for k = 1:numel(cpi_span)    
    EFTDist(3*k-2) = DOFNumbering(cpi_span(k),cpj_span(k),1);
    EFTDist(3*k-1) = DOFNumbering(cpi_span(k),cpj_span(k),2);
    EFTDist(3*k)   = DOFNumbering(cpi_span(k),cpj_span(k),3);        
end 

for elj = q+1:meta-q-1
    for eli = p+1:mxi-p-1
        % check if element is greater than zero
        if Xi(eli+1)~=Xi(eli) && Eta(elj+1)~=Eta(elj)
            %% 3i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
            %
            %         | xi_i+1 - xi_i                    |
            %         | -------------            0       |
            %         |        2                         |
            %  xi,u = |                                  |
            %         |                  eta_j+1 - eta_j |
            %         |        0         --------------- |
            %         |                          2       |
            detJxiu = (Xi(eli+1)-Xi(eli))*(Eta(elj+1)-Eta(elj))/4;
            
            %% 3ii. Create the Element Freedom Table
            
            % Initialize element freedom table
            EFT = zeros(1,noDOFsEl);
            % initialize counter
            k = 1;
                                     
            % relation global-local dof
            for cpj = elj-q:elj
                for cpi = eli-p:eli
                    EFT(k) = DOFNumbering(cpi,cpj,1);
                    EFT(k+1) = DOFNumbering(cpi,cpj,2);
                    EFT(k+2) = DOFNumbering(cpi,cpj,3);
                    % Update counter
                    k = k + 3;                                    
                end
            end
            
            % check if the element holds the disturbed DOF (i,j)
            if any((eli-p:eli)==disturbed_cp(1)) && ...
                any((elj-q:elj)==disturbed_cp(2))
                flag=true;
            else
                flag=false;
            end
                      
            %% 3iii. Initialize the element area
            elementArea = 0;
            elementAreaDist=0;
            %% 3iv. Loop over all the Gauss Points
            for cEta = 1:etaNGP
                for cXi = 1:xiNGP
                    %% 3iv.1. Compute the NURBS coordinates u,v of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                    xi = ( Xi(eli+1)+Xi(eli) + xiGP(cXi)*(Xi(eli+1)-Xi(eli)) )/2;
                    eta = ( Eta(elj+1)+Eta(elj) + etaGP(cEta)*(Eta(elj+1)-Eta(elj)) )/2;
                    
                    %% 3iv.2. Compute the NURBS basis function and up to their second derivatives at the Gauss Point
                    nDrvBasis = 2;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface(eli,p,xi,Xi,elj,q,eta,Eta,CP,isNURBS,nDrvBasis);                    
                    
                    %% 3iv.3. Compute the covariant base vectors and their first derivatives
                    nDrvBaseVct = 1;
                    [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface(eli,p,elj,q,CP,nDrvBaseVct,dR);
                 
                    %% 3iv.4. Compute the surface normal (third covariant base vector not normalized)
                    G3Tilde = cross(dG1(:,1),dG2(:,1));
                    
                    %% 3iv.5. Compute the legth of G3Tilde (= area dA)
                    dA = norm(G3Tilde);
  
                    %% 3iv.6. Compute the element stiffness matrix at the Gauss point
                    KeOnGP = computeElStiffMtxKirchhoffLoveShellLinear(p,q,dR,[dG1(:,1) dG2(:,1)],[dG1(:,2) dG2(:,2) dG1(:,3)],G3Tilde,Dm,Db);
                    
                    %% 3iv.7 Compute the element area on the Gauss Point and add the contribution
                    elementAreaOnGP = dA*detJxiu*xiGW(cXi)*etaGW(cEta);
                    elementArea = elementArea + elementAreaOnGP;                    
                   
                    %% 3iv.8. Add the contribution from the Gauss Point
                    GaussContribution=KeOnGP*elementAreaOnGP;
                    K(EFT,EFT) = K(EFT,EFT) + GaussContribution;
                                        
                    %% Compute distorted Stiffness matrix only if the flag is active
                    if flag==true % recompute stiffness matrix with distorted element
                        dRDist = computeIGABasisFunctionsAndDerivativesForSurface(eli,p,xi,Xi,elj,q,eta,Eta,CPDist,isNURBS,nDrvBasis);
                        [dG1Dist,dG2Dist] = computeBaseVectorsAndDerivativesForBSplineSurface(eli,p,elj,q,CPDist,nDrvBaseVct,dRDist);
                        G3TildeDist = cross(dG1Dist(:,1),dG2Dist(:,1));
                        dADist = norm(G3TildeDist);
                        KeOnGPDist = computeElStiffMtxKirchhoffLoveShellLinear(p,q,dRDist,[dG1Dist(:,1) dG2Dist(:,1)],[dG1Dist(:,2) dG2Dist(:,2) dG1Dist(:,3)],G3TildeDist,Dm,Db);
                        elementAreaOnGPDist = dADist*detJxiu*xiGW(cXi)*etaGW(cEta);
                        elementAreaDist = elementAreaDist + elementAreaOnGPDist;
                        %computes the stiffness matrix
                        KDist(EFT,EFT) = KDist(EFT,EFT) + KeOnGPDist*elementAreaOnGPDist;
                    else
                        KDist(EFT,EFT) = KDist(EFT,EFT) + GaussContribution;
                    end
                end
            end
            %% 3v. Find the minimum element area in the mesh
            if elementArea<minElArea
                minElArea = elementArea;
            end
        end
    end
end

F = zeros(noDOFs,1);

%% 4. Compute the exernally applied load vector
t = 0;    
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    F = funcHandle(F,NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},p,q,Xi,Eta,CP,isNURBS,...
        NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection(counterNBC,1),t,int,'');
end

end