function [ Ep_final, delta_final, RelErr, Ep_history, delta_history,MassFinal ] = SensitivityWithErrorChecks( BSplinePatch,K,CP2Dist,vector,dHatLinear,RelErrTolerance,delta_in )
%SENSITIVITYWITHERRORCHECKS Calculates the sensitivity for a given

%disturbance in the control points of a BSplinePatch.
%Cases are adjusted for counting as well the new parameters of total mass
%and min element area
switch nargin
    case 5
        RelErrTolerance = 10^-5;
        delta_in = -1; 
    case 6 
        delta_in = -1; 
    case 7
        %do nothing
    otherwise
        error('not enough input arguments!');
end
        
%% Solve the system applying linear analysis

RelErr = inf; % initial relative error
EpOld = inf; % initial EpOld
MassDist=0; % mass of the disturbed model
if delta_in < 0 % use default delta and find optimal delta iteratively
    delta = sqrt(BSplinePatch.minElArea); % initial delta for finite differences (scaled to min element size)
    do_not_iterate = 0; % iterate until error bounds are met
else % use given delta
    delta = delta_in; % use input delta
    do_not_iterate = 1; % just do one iteration, then break
end

Ep_history = [];
iteration_count = 0;  

while (RelErr > RelErrTolerance) % check if error meets the demands of the user or iteration is prohibited
    iteration_count = iteration_count + 1;
    
    %returns the BSPLINEPATCH with the modified control points stored in the
    %variable CPd
    
    [BSplinePatch]=CPDisturbance(BSplinePatch,CP2Dist,vector,delta,1);      
    
    [KDist, dindex,MassDist] = computeLinearMtrcsSensitivity(BSplinePatch,CP2Dist);        
    
    [ Ep ] = Sensitivity(K,KDist,delta,dHatLinear,dindex);
     
    Ep_history(iteration_count)=Ep;
    delta_history(iteration_count)=delta;
    RelErr = abs((Ep - EpOld)/Ep);
    EpOld = Ep;
    delta=delta/2;%this is the increment used to calculate the new CP, value that converged for several tests
    if (delta < 10^6*eps) % output error message if delta goes too low, below machine prescision!
        warning(['sensitivity analysis has not converged up to the given relative error tolerance of ',mat2str(RelErrTolerance),'!\n current error: ',mat2str(RelErr)]);
        break;
    end
    if do_not_iterate % if iteration is prohibited just take the first iteration for the final result
        break;
    end
end

Ep_final = Ep_history(end);
MassFinal= (MassDist-BSplinePatch.TotalMass) / delta;%calculates the sensitivity directly
delta_final = delta_history(end);

end

