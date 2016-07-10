function [ Ep_final, Mass_final, Disp_final, RelErr, delta ] = SensitivityWithErrorChecks( BSplinePatch,CP2Dist,vector,K_global,u_global,Compute,RelErrTolerance,delta_in )
%SENSITIVITYWITHERRORCHECKS Calculates the sensitivity for a given
Mass_final=0;
Ep_final=0;
Disp_final=0;
%disturbance in the control points of a BSplinePatch.
%Cases are adjusted for counting as well the new parameters of total mass
%and min element area

if nargin < 3
    error('not enough input arguments!');
end

if nargin == 4
    error('please provide K_global and u_global.');
end

if nargin < 5
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
    %Precompute the element stiffness matrix of the unperturbed model and save
    %them to the BSplinePatch
    [BSplinePatch] = computeElementStiffnessMatrices(BSplinePatch);
    [K_global, F_global] = assembleGlobalSystem(BSplinePatch);
    u_global = solveGlobalSystem(K_global, F_global, BSplinePatch, solve_LinearSystem);
end

if nargin < 6
    Compute = 'all';
end

if nargin < 7    
    RelErrTolerance = 10^-5;
end

if nargin < 8
   delta_in = -1; 
end
        
%% checks which sensitivity will be computed

if strcmp(Compute,'displacement')
    disp_flag=true;
    en_flag=false;
    mass_flag=true;
elseif strcmp(Compute,'energy')
    disp_flag=false;
    en_flag=true;
    mass_flag=true;
elseif strcmp(Compute,'mass')
    disp_flag=false;
    en_flag=false;
    mass_flag=true;
elseif strcmp(Compute,'all')
    disp_flag=true;
    en_flag=true;
    mass_flag=true;
else    
    error('unknown key for Compute. Please provide one of the following: displacement, energy, mass or all.');
end        
    
%% Solve the system applying linear analysis

RelErr = inf; % initial relative error
EpOld = inf; % initial EpOld
EpDispOld = inf; % initial EpDispOld
MassOld = inf; %initial MassDistOld

if delta_in < 0 % use default delta and find optimal delta iteratively
    delta = sqrt(BSplinePatch.minElArea)*2; % initial delta for finite differences (scaled to min element size)
    do_not_iterate = 0; % iterate until error bounds are met
else % use given delta
    delta = delta_in*2; % use input delta
    do_not_iterate = 1; % just do one iteration, then break
end

Ep_history = [];
iteration_count = 0;  

while (RelErr > RelErrTolerance) % check if error meets the demands of the user or iteration is prohibited
    delta=delta/2;%this is the increment used to calculate the new CP, value that converged for several tests
    
    iteration_count = iteration_count + 1;
    
    %returns the BSPLINEPATCH with the modified control points stored in the
    %variable CPd    
    [BSplinePatch] = CPDisturbance(BSplinePatch,CP2Dist,vector,delta,1);      
    
    [KDist, dindex, MassDist] = computeLinearMtrcsSensitivity(BSplinePatch,CP2Dist);                
    
    RelErr_Ep = 0;
    RelErr_EpDisp = 0;
    RelErr_Mass = 0;
    
    %figure(1)
    %hold on
    
    if en_flag == true
        Ep = Sensitivity(K_global,KDist,delta,u_global,dindex);        
        RelErr_Ep = abs((Ep - EpOld)/Ep);
        EpOld = Ep;
        %plot(iteration_count,RelErr_Ep,'r*');
    end
    if disp_flag == true
        EpDisp = DisplacementSensitivity(K_global,KDist,delta,u_global,dindex,BSplinePatch);        
        RelErr_EpDisp = abs((EpDisp - EpDispOld)/EpDisp);
        EpDispOld = EpDisp;
        %plot(iteration_count,RelErr_EpDisp,'b*');
    end
    if mass_flag == true
        Mass = (MassDist-BSplinePatch.TotalMass) / delta;
        RelErr_Mass = abs((Mass - MassOld)/Mass);
        MassOld = Mass;
        %plot(iteration_count,RelErr_Mass,'k*');
    end   
    
    %set(gca,'YScale','log');
    
    %hold off
    
    RelErr = max([RelErr_Ep, RelErr_EpDisp, RelErr_Mass]);
        
    if (delta < 10^6*eps) % output error message if delta goes too low, below machine prescision!
        warning(['sensitivity analysis has not converged up to the given relative error tolerance of ',mat2str(RelErrTolerance),'!\n current error: ',mat2str(RelErr)]);
        break;
    end
    if do_not_iterate % if iteration is prohibited just take the first iteration for the final result
        break;
    end
end
if en_flag
    Ep_final = Ep;
end
if disp_flag
    Disp_final = EpDisp;
end
if mass_flag==true
    Mass_final= Mass;
end

end

