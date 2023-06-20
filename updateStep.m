 function [LSgridPhi2] = updateStep(LSgridPhi2,shapeSens,stepLength,iterNum)
 % Smooth the sensitivities
 %[shapeSens] = conv2(padarray(shapeSens,[1,1],'replicate'),1/6*[0 1 0; 1 2 1;0 1 0],'valid');
 % Load bearing pixels must remain solid -Bridge:
 
 % Design update via evolution
 [LSgridPhi2] = evolve(-shapeSens,LSgridPhi2,stepLength);