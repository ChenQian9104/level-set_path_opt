function [SignDistPhi] = Reinitialize(Phi0, dx, dy, LoopNum)
%=================================================================%
%=================================================================%
for i = 1 : LoopNum + 1
    [Dx_L, Dx_R] = UpwindDiff(Phi0 , dx , 'x');
    [Dy_L, Dy_R] = UpwindDiff(Phi0 , dy , 'y');
    Dx_C = (Dx_L + Dx_R)/2;
    Dy_C = (Dy_L + Dy_R)/2;
    S = Phi0 ./ (sqrt(Phi0.^2 + (Dx_C.^2 + Dy_C.^2) * dx^2) + eps);
    DetT = 0.5 * min(dx,dy)/max(abs(S(:)));
    Grad_Plus  = ((max(Dx_L,0)).^2 + (min(Dx_R , 0)).^2 + (max(Dy_L,0)).^2 + (min(Dy_R,0)).^2 ).^0.5;
    Grad_Minus = ((min(Dx_L,0)).^2 + (max(Dx_R , 0)).^2 + (min(Dy_L,0)).^2 + (max(Dy_R,0)).^2 ).^0.5;
    Phi0 = Phi0 - DetT .* ((max(S, 0) .* Grad_Plus + min(S, 0) .* Grad_Minus) - S);
end;
SignDistPhi = Phi0;