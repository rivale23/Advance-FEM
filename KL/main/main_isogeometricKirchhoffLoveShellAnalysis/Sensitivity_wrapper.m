function [ SensitivityMatrix ] = Sensitivity_wrapper( BSplinePatch, vectors )
%SENSITIVITY_WRAPPER Summary of this function goes here
%   Detailed explanation goes here

RelErrTolerance = 10^(-5);
SensitivityMatrix = zeros([size(vectors),3]);

solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% compute solution once for all sensitivities
[dHatLinear,~,~,~] = solve_IGAKirchhoffLoveShellLinear_shortcut...
    (BSplinePatch,[],solve_LinearSystem);

delta = -1; % initial finite difference delta equal to -1 allows iteration for the first CP, other control points use this delta to save computation time

for i = 1:size(vectors,1)
    for j = 1:size(vectors,2)
        vector=vectors{i,j};%direction of the distortion
        for d = 1:3
            vector_component = zeros(size(vector));
            vector_component(d) = vector(d);
            if max(abs(vector_component)) == 0
                disp(['sensitivity for component ',mat2str(d),' of CP @',mat2str([i,j]),' not computed, because disturbance equal to zero.']);
                SensitivityMatrix(i,j,d)=0;
            else
                disp(['calculating sensitivity for component ',mat2str(d),' of CP @ ',mat2str([i,j])]);                
                CP2Dist=[i j];%control pint to disturb         
                [Ep_final, delta, ~, ~, ~] = SensitivityWithErrorChecks( BSplinePatch,CP2Dist,vector_component,dHatLinear,RelErrTolerance,delta);                
                SensitivityMatrix(i,j,d) = Ep_final; % save sensitivity to matrix
            end
        end
    end
end

end

