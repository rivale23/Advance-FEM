function [ SensitivityMatrix ] = SensitivityWrapper2( BSplinePatch, vectors, DesiredCPs,flag )
%SENSITIVITY_WRAPPER Summary of this function goes here
%If the flag is set to 1, then the sensitivity will be computed independen
%for each direction X,Y and Z, otherwise, the whole vector will be applied
%as disturbance with a value of delta

%if the flag is not given
if nargin < 4
    flag=false;
end
sv=size(vectors,2);
RelErrTolerance = 10^(-5);
%the falg will 
if flag==true
    SensitivityMatrix=zeros([sv,3]);  
else
    SensitivityMatrix=zeros([sv,1]); 
end
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
BSplinePatch = computeElementStiffnessMatrices(BSplinePatch);
[K_global, F_global] = assembleGlobalSystem(BSplinePatch);
dHatLinear = solveGlobalSystem(K_global, F_global, BSplinePatch, solve_LinearSystem);

delta = -1; % initial finite difference delta equal to -1 allows iteration for the first CP, other control points use this delta to save computation time

for i = 1:sv
    vector=vectors{i};
    CP2Dist=DesiredCPs{i};
    if flag==true
        for d = 1:3
            vector_component = zeros(size(vector));
            vector_component(d) = vector(d);
            if max(abs(vector_component)) == 0
                disp(['sensitivity for component ',mat2str(d),' of CP @',mat2str([i,j]),' not computed, because disturbance equal to zero.']);
                SensitivityMatrix(i,d)=0;
            else    
                [Ep_final, delta, a, b, c] = SensitivityWithErrorChecks( BSplinePatch,K_global,CP2Dist,vector_component,dHatLinear,RelErrTolerance,delta);                
                SensitivityMatrix(i,d)=Ep_final;
            end
        end
    else
        [Ep_final, delta, a, b, c] = SensitivityWithErrorChecks( BSplinePatch,K_global,CP2Dist,vector,dHatLinear,RelErrTolerance,delta);                
        SensitivityMatrix(i)=Ep_final;
    end
end
end