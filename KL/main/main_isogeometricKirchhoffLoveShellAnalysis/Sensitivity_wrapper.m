function [ SensitivityMatrix ] = Sensitivity_wrapper( BSplinePatch, vectors )
%SENSITIVITY_WRAPPER Summary of this function goes here
%   Detailed explanation goes here

RelErrTolerance = 10^(-5);
SensitivityMatrix = zeros(size(vectors));

solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

for i = 1:size(vectors,1)
    for j = 1:size(vectors,2)
        disp(['calculating sensitivity for CP @ ',mat2str([i,j])]);
        vector=vectors{i,j};%direction of the distortion
        CP2Dist=[i j];%control pint to disturb         
        [Ep_final, ~, ~, ~, ~] = SensitivityWithErrorChecks( BSplinePatch,CP2Dist,vector,solve_LinearSystem,RelErrTolerance );
        SensitivityMatrix(i,j) = Ep_final; % save sensitivity to matrix
    end
end    

end

