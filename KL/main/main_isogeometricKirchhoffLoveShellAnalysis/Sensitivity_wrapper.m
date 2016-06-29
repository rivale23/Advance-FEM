function [ SensitivityMatrix,SensitivityMass ] = Sensitivity_wrapper( BSplinePatch, vectors, IndependentDirectionsFlag)
%SENSITIVITY_WRAPPER Summary of this function goes here
%if IndependentDirectionsFlag is send as true, then it will compute the
%sensitivity indepenedent for each directions (X,Y,Z), if the argument is false or
%is not given, it will compute only one sensitivity (direction of the vector)
%it returns the sentitivities wrt to the strain energy and the mass
%function for the requested control points and directions


%% check input and add default values if input is missing

if nargin < 3
    IndependentDirectionsFlag=false;
end

if IndependentDirectionsFlag==true
   SensitivityMatrix = zeros([size(vectors),3]);
   SensitivityMass = zeros([size(vectors),3]);
else
    SensitivityMatrix = zeros([size(vectors),1]);
    SensitivityMass = zeros([size(vectors),1]);
end

%% precompute static variables like displacement and stiffness matrices

RelErrTolerance = 10^(-5);

solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
%Precompute the element stiffness matrix of the unperturbed model and save
%them to the BSplinePatch
[BSplinePatch] = computeElementStiffnessMatrices(BSplinePatch);
[K_global, F_global] = assembleGlobalSystem(BSplinePatch);
u_global = solveGlobalSystem(K_global, F_global, BSplinePatch, solve_LinearSystem);

delta = -1; % initial finite difference delta equal to -1 allows iteration for the first CP, other control points use this delta to save computation time

%% iterate over all displacement directions and compute sensitivities
for i = 1:size(vectors,1)
    for j = 1:size(vectors,2)        
        vector=vectors{i,j};%direction of the distortion
        CP2Dist=[i j];%control point to disturb
        if IndependentDirectionsFlag==true % compute all sensitivities for each vector component independently
            for d = 1:3
                vector_component = zeros(size(vector));
                vector_component(d) = vector(d);
                if max(abs(vector_component)) == 0
                    disp(['sensitivity for component ',mat2str(d),' of CP @',mat2str([i,j]),' not computed, because disturbance equal to zero.']);
                    SensitivityMatrix(i,j,d)=0;
                    SensitivityMass(i,j,d)=0;
                else
                    disp(['calculating sensitivity for component ',mat2str(d),' of CP @ ',mat2str([i,j])]);  
                    %calculation of the sensitivity of the strain energy
                    %and mass. The minimum size area is passed for the
                    %initial guess of the default convergence analysis of
                    %the derivatives, which is done only for the first
                    %requested control point, for efficiency reasons. It
                    %also returns the sensitivity of mass of the model due
                    %to the perturbation
                    [Ep_final, delta, ~, ~, ~, Mass_final] = SensitivityWithErrorChecks( BSplinePatch,CP2Dist,vector_component,K_global,u_global,RelErrTolerance,delta);                
                    SensitivityMatrix(i,j,d) = Ep_final; % save sensitivity to matrix
                    %The mass is currently calculated assuming a constant density of
                    %1.0 and constant thickness. Here the thickness is
                    %multiplied to the total mass (area) of the shell. By
                    %the definition of derivative,sensitivity=(final_mass-initial_mass)/delta = 
                    %= 1.0*thickness*(final_area-initial_area)delta , then 
                    %it results the same to multiply the thickness to the sensitivity 
                    %than to the perturbed and unperturbed shells
                    SensitivityMass(i,j,d) = Mass_final;%save sensitivity mass
                end
            end
        else % only compute sensitivity for each vector one with all three components combined
            if max(abs(vector)) == 0
                disp(['sensitivity for CP @',mat2str([i,j]),' not computed, because disturbance equal to zero.']);
                SensitivityMatrix(i,j)=0;
                SensitivityMass(i,j)=0;
            else
                disp(['calculating sensitivity of CP @ ',mat2str([i,j])]);
                [Ep_final, delta, ~, ~, ~, Mass_final] = SensitivityWithErrorChecks( BSplinePatch,CP2Dist,vector,K_global,u_global,RelErrTolerance,delta);                
                SensitivityMatrix(i,j) = Ep_final; % save sensitivity to matrix
                SensitivityMass(i,j) = Mass_final;%save sensitivity mass
            end
        end
    end
end

end

