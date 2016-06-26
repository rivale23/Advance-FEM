function [KDist, EFTDist, MassDist] = computeLinearMtrcsSensitivity(BSplinePatch,disturbed_cp)


%% 0. Read input

if ~isfield(BSplinePatch,'elements') % element stiffness matrices have not been precomputed.
    BSplinePatch = computeElementStiffnessMatrices(BSplinePatch);
end

% Assign back the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CPDist= BSplinePatch.CPd;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
int = BSplinePatch.int;
DOFNumbering = BSplinePatch.DOFNumbering;
MassDist=0;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
% Number of control points in xi-,eta-direction
nxi = length(CPDist(:,1,1));
neta = length(CPDist(1,:,1));

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Initialize minimum element area in the IGA mesh
tolerance = 1e-4;
if abs(CPDist(1,1,1)-CPDist(nxi,1,1))>=tolerance
    minElArea = abs(CPDist(1,1,1)-CPDist(nxi,1,1));
else
    minElArea = CPDist(1,1,1)-CPDist(1,neta,1);
end

% Global number of DOFs
noDOFs = 3*nxi*neta;

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
        %% 3i. Read precomputed quantities        
        element = BSplinePatch.elements{eli-p,elj-q};                                
        % Read element freedom table
        EFT = element.EFT;                                       
        
        %% 3ii. check if the element holds the disturbed DOF (i,j)
        if any((eli-p:eli)==disturbed_cp(1)) && ...
            any((elj-q:elj)==disturbed_cp(2))
            flag=true;
        else
            flag=false;
        end 
        
        %% 3iii. reuse or recompute element stiffness matrix                
        if flag % only recompute element stiffness matrix, if distortion on element occours                
            [ K_local, elementAreaDist ] = computeElementStiffnessMatrix(eli,p,Xi,elj,q,Eta,CPDist,isNURBS,xiNGP,xiGP,xiGW,etaNGP,etaGP,etaGW,Dm,Db);                                            
        else % otherwise reuse precomputed element matrix
            elementAreaDist = element.elementArea;
            K_local = element.K_local;                
        end
        MassDist=MassDist+elementAreaDist;
        KDist(EFT,EFT) = KDist(EFT,EFT) + K_local;

        %% 3iv. Find the minimum element area in the mesh
        if elementAreaDist<minElArea
            minElArea = elementAreaDist;
        end        
    end
end

MassDist = MassDist * BSplinePatch.parameters.t;

end