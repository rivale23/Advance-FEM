function [ Es ] = Sensitivity( Kdelta,delta,U )
%SENSITIVITY computes the sensitivity of the desired Control points
%gets the array of Kdelta and for the several control points with the whole
%dimension of the matrix and gets the whole U. and will return an array of
%all the sensitivities Es

%reduces each matrix by deleting all the zero entries
    Kdim=any(Kdelta,1);
    %get the indices of non zero elements
    Kindex=find(Kdim);
    %resize the matrix
    KD=Kdelta(Kindex,Kindex);
    %computes the sensitivity
    
    Ured=U(Kindex);
    Kp=(KD)./delta;

    Es=((-Kp*Ured)'*Ured)./2;

end

