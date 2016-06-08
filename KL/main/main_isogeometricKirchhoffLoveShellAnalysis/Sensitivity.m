function [ Ep ] = Sensitivity( Kdelta,delta,U )
%SENSITIVITY Summary of this function goes here
%   Detailed explanation goes here

Kp=(Kdelta)./delta;

Ep=((-Kp*U)'*U)./2;

end

