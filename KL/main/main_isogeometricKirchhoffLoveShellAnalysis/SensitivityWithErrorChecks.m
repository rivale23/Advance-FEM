function [ Ep_final, RelErr, Ep_history, delta_history, dHatLinear ] = SensitivityWithErrorChecks( BSplinePatch,CP2Dist,vector,solve_LinearSystem,RelErrTolerance )
%SENSITIVITYWITHERRORCHECKS Calculates the sensitivity for a given
%disturbance in the control points of a BSplinePatch.

%% %Advance FEM:  Solve the system applying linear analysis

RelErr = inf; % initial relative error
EpOld = inf; % initial EpOld
delta = 0.01; % initial delta for finite differences
iteration_count = 0;

Ep_history = [];

while (RelErr > RelErrTolerance) % check if error meets the demands of the user
    iteration_count = iteration_count + 1;
    
    %returns the BSPLINEPATCH with the modified control points stored in the
    %variable CPd
    [BSplinePatch]=CPDisturbance(BSplinePatch,CP2Dist,vector,delta,1);
    [KDist,K,dindex]=ReducedStiffnessMatrix(BSplinePatch,CP2Dist);

    %% Attention!!!! the function ReducedStiffnessMatrix needs to replace the original function, send all the requiered arguments
    %% -------------------------------------
    [dHatLinear,F,minElArea,StiffnessMatrix] = solve_IGAKirchhoffLoveShellLinear...
        (BSplinePatch,solve_LinearSystem,'');

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
end

Ep_final = Ep_history(end);


end

