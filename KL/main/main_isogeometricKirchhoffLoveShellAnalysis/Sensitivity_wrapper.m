function [ SensitivityMatrix ] = Sensitivity_wrapper( BSplinePatch, vectors, IndependentDirectionsFlag )
%SENSITIVITY_WRAPPER Summary of this function goes here
%if IndependentDirectionsFlag is send as true, then it will compute the
%sensitivity indepenedent for each directions (X,Y,Z), if the argument is false or
%is not given, it will compute only one sensitivity (direction of the vector)

if nargin<3
    IndependentDirectionsFlag=false;
end

if IndependentDirectionsFlag==true
   SensitivityMatrix = zeros([size(vectors),3]);
else
    SensitivityMatrix = zeros([size(vectors),1]);
end
RelErrTolerance = 10^(-5);

solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

BSplinePatch = computeElementStiffnessMatrices(BSplinePatch);
[K_global, F_global] = assembleGlobalSystem(BSplinePatch);
dHatLinear = solveGlobalSystem(K_global, F_global, BSplinePatch, solve_LinearSystem);

delta = -1; % initial finite difference delta equal to -1 allows iteration for the first CP, other control points use this delta to save computation time

for i = 1:size(vectors,1)
    for j = 1:size(vectors,2)        
        vector=vectors{i,j};%direction of the distortion
        CP2Dist=[i j];%control pint to disturb
        if IndependentDirectionsFlag==true
            for d = 1:3
                vector_component = zeros(size(vector));
                vector_component(d) = vector(d);
                if max(abs(vector_component)) == 0
                    disp(['sensitivity for component ',mat2str(d),' of CP @',mat2str([i,j]),' not computed, because disturbance equal to zero.']);
                    SensitivityMatrix(i,j,d)=0;
                else
                    disp(['calculating sensitivity for component ',mat2str(d),' of CP @ ',mat2str([i,j])]);                                        
                    [Ep_final, delta, a, b, c] = SensitivityWithErrorChecks( BSplinePatch,K_global,CP2Dist,vector_component,dHatLinear,RelErrTolerance,delta);                
                    SensitivityMatrix(i,j,d) = Ep_final; % save sensitivity to matrix
                end
            end
        else
            disp(['calculating sensitivity of CP @ ',mat2str([i,j])]);
            [Ep_final, delta, a, b, c] = SensitivityWithErrorChecks( BSplinePatch,K_global,CP2Dist,vector,dHatLinear,RelErrTolerance,delta);                
            SensitivityMatrix(i,j) = Ep_final; % save sensitivity to matrix
        end
    end
end

end

