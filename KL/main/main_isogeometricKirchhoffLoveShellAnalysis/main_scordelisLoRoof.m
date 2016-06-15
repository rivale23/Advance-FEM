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
%% Script documentation
% 
% Task : The full Scordelis-Lo-Roof benchmark problem is modelled and
%        solved in a single patch geometry.
%
% Date : 12.02.2015
%
%% Preamble
clear; 
clc;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Include linear equation system solvers
addpath('../../equationSystemSolvers/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric Kirchhoff-Love shell formulation
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/');

%% NURBS parameters

% Global variables
Length = 50;
Radius = 25;

% Polynomial degrees
p = 1;
q = 2;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-Length/2 -Length/2 -Length/2
             Length/2  Length/2  Length/2];
         
% y-coordinates
CP(:,:,2) = [-Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)
             -Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)];
         
% z-coordinates
CP(:,:,3) = [Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)
             Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)];
       
% Weights
weight = cos(2*pi/9);
CP(:,:,4) = [1 weight 1
             1 weight 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j=1:neta
        if CP(i,j,4)~=1
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants

% Young's modulus
parameters.E = 4.32e8;

% Poisson ratio
parameters.nue = .0;

% Thickness of the shell
parameters.t = .25;

% Density of the shell (used only for dynamics)
parameters.rho = 7850;

%% GUI

% Analysis type
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'manual')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
end

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment'
graph.resultant = 'force';

% Component of the resultant to plot
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2Principal';

%% Refinement

% Degree by which to elevate
tp = 1;
tq = 0;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 1;
edgeRatio = ceil(Length/Radius/(sin(4*pi/9)));
refXi = edgeRatio*scaling;
refEta = ceil(4/3)*scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)

% back and front curved edges are a rigid diaphragm
homDOFs = [];
xiSup = [0 0];   etaSup = [0 1];    
for dirSupp = [2 3];
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];    
for dirSupp = [2 3];
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix the back left corner of the shell to avoid rigid body motions
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 1;
homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC = [];

% Cables
cables.No = 0;

% load (Neuman boundary conditions)
FAmp = - 9e1;
NBC.noCnd = 1;
xib = [0 1];   etab = [0 1];   dirForce = 3;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection(1,1) = dirForce;
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isConservative(1,1) = true;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,valuesInhomDOFs,...
    weakDBC,cables,NBC,[],[],[],[],[],int);


%% Advance FEM: here The control points can be modified to change the shape here.
magnitude=0
vector=[1,0,0];%direction of the distortion
CP2Dist=[5 1];%control pint to disturb
[BSplinePatch]=CPDisturbance(BSplinePatch,CP2Dist,vector,magnitude,0);


%% Compute the load vectors for each patch (only for the visualization)
FGamma = zeros(3*BSplinePatch.noCPs,1);
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    FGamma = funcHandle...
        (FGamma,NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},BSplinePatch.p,...
        BSplinePatch.q,BSplinePatch.Xi,BSplinePatch.Eta,...
        BSplinePatch.CP,BSplinePatch.isNURBS,NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection(counterNBC,1),0,BSplinePatch.int,'outputEnabled');
end
    


%% Plot reference configuration
figure(graph.index)
%Advance FEM: to plot the reference shell I just change the argument from
%CP to BSplinePatch.CP
plot_referenceConfigurationIGAThinStructure(p,q,Xi,Eta,BSplinePatch.CP,isNURBS,homDOFs,FGamma,'outputEnabled');
title('Reference configuration for an isogeometric Kirchhoff-Love shell');
graph.index = graph.index + 1;

%% %Advance FEM:  Solve the system applying linear analysis
%------------------------------------
%sensitivity analysis
%modify the control points as desired

delta = 0.01;
iteration_count = 0;
RelErr = inf;
RelErrTolerance = 10^(-5);
EpOld = inf;

while (RelErr > RelErrTolerance)
    iteration_count = iteration_count + 1
    
    %returns the BSPLINEPATCH with the modified control points stored in the
    %variable CPd
    [BSplinePatch]=CPDisturbance(BSplinePatch,CP2Dist,vector,delta,1);


    [KDist,K,dindex]=ReducedStiffnessMatrix(BSplinePatch,CP2Dist);




    %% Attention!!!! the function ReducedStiffnessMatrix needs to replacethe original function, send all the requiered arguments
    %% -------------------------------------
    [dHatLinear,F,minElArea,StiffnessMatrix] = solve_IGAKirchhoffLoveShellLinear...
        (BSplinePatch,solve_LinearSystem,'');

    [ Ep ] = Sensitivity(K,KDist,delta,dHatLinear,dindex);
    EpV(iteration_count)=Ep;
    deltaV(iteration_count)=delta;
    RelErr = abs((Ep - EpOld)/Ep);
    EpOld = Ep;
    delta=delta/2;%this is the increment used to calculate the new CP, value that converged for several tests
    if (delta < 10^6*eps)
        warning(['sensitivity analysis has not converged up to the given relative error tolerance of ',mat2str(RelErrTolerance),'!\n current error: ',mat2str(RelErr)]);
        break;
    end
end
figure(9)
hold on
semilogx(1./deltaV,EpV);
xlabel('1/delta');
ylabel('Sensitivity');
hold off

figure(12)
hold on
t3=KDist-K;
surface(t3);
hold off

%% Postprocessing
graph.index = plot_postprocIGAKirchhoffLoveShellLinear(BSplinePatch,dHatLinear,graph,'outputEnabled');
title('Linear analysis');
