function [BackDiff, FawdDiff] = UpwindDiff(Phi , dx , strDirection)
%=================================================================%
% function [BkDiff, FwDiff] = UpwindDiff(Phi , dx) calculates
% backward and forward finite difference,where
% Phi: is an n-by-n matrix;
% dx : the interval between two adjacent grids in axis X.
% strDirection : is a character string. It equals to 'x'or 'y'which mean
% get spacial derivatives in x direction or y direction.
% BackDiff:is an n-by-n matrix, which stores (Phi(i,j) - Phi(i-1 ,j))/dx;
% FawdDiff:is an n-by-n matrix, which stores (Phi(i+1,j) - Phi(i ,j))/dx;
%   
%=================================================================%
Matrix = Matrix4diff( Phi );
if strDirection == 'x'
    BackDiff = (Phi - Matrix.i_minus_1)/dx;
    FawdDiff = (Matrix.i_plus_1 - Phi)/dx;
elseif strDirection == 'y'
    BackDiff = (Phi - Matrix.j_minus_1)/dx;
    FawdDiff = (Matrix.j_plus_1 - Phi)/dx;  
end;