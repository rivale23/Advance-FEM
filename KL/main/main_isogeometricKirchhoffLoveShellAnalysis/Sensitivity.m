function [ Ep ] = Sensitivity( K,Kdelta,delta,U )
%SENSITIVITY Summary of this function goes here
%   Detailed explanation goes here

Kp=(Kdelta-K)./delta;

Ep=((-Kp*U)'*U)./2;

end

