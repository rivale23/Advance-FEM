function [ Es ] = DisplacementSensitivity( Ko,Kdelta,delta,U, Kindex,BSplinePatch)
%SENSITIVITY computes the sensitivity of the desired Control points
%gets the array of Kdelta and for the several control points with the whole
%dimension of the matrix and gets the whole U. and will return an array of
%all the sensitivities Es

solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
%reduces each matrix by deleting all the zero entries
    %resize the matrices
    KD=Kdelta;
    K=Ko;
    %computes the sensitivity
    
    Ured=U;
    Kp=(KD-K)./delta;
    V=zeros(size(U));
    V(Kindex)=V(Kindex)+1;
    rhs = Kp*Ured;
    
    %needs to solve the sytem considering the boundary conditions, otherwise
    %the matrix is singular
    lhs = solveGlobalSystem(K, rhs, BSplinePatch, solve_LinearSystem);
    Es=-(lhs)' ;   
    Es=Es*V;    
end 
