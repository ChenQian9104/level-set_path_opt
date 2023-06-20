function [Matrix] = Matrix4diff( Phi )
%=================================================================%
% function [Matrix] = Matrix4diff( Phi ) produces a structure used for
% upwind finite diffence.
% Phi: is an m-by-n matrix;
% Matrix: a structure which includes 4 matrixes used to calculate finite
% difference.
%=================================================================%
Matrix.i_minus_1 = zeros(size(Phi));
Matrix.i_plus_1 = Matrix.i_minus_1; 
Matrix.j_minus_1 = Matrix.i_minus_1; 
Matrix.j_plus_1 = Matrix.i_minus_1;
Matrix.i_minus_1(:, 1) = Phi(:, end); 
Matrix.i_minus_1(:, 2:end) = Phi(:,1:end-1);
Matrix.i_plus_1(:, end) = Phi(:,1); 
Matrix.i_plus_1(:, 1:end-1) = Phi(:,2:end);
Matrix.j_minus_1(1, :) = Phi(end, :); 
Matrix.j_minus_1(2:end , :) = Phi(1:end-1,:);
Matrix.j_plus_1(end,:) = Phi(1,:); 
Matrix.j_plus_1(1:end-1, :) = Phi(2:end, :);