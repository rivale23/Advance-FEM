function [ Es ] = DisplacementSensitivity( Ko,Kdelta,delta,U, Kindex)
%SENSITIVITY computes the sensitivity of the desired Control points
%gets the array of Kdelta and for the several control points with the whole
%dimension of the matrix and gets the whole U. and will return an array of
%all the sensitivities Es

%reduces each matrix by deleting all the zero entries
    %resize the matrices
    KD=Kdelta(Kindex,Kindex);
    K=Ko(Kindex,Kindex);
    %computes the sensitivity
    
    Ured=U(Kindex);
    Kp=(KD-K)./delta;
    V=Kindex;
    rhs = Kp*Ured;
    Es = -(K\rhs)';    
    Es=Es*V'./2;    
end 
