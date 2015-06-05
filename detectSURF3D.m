function [points] = detectSURF3D(V)
% DETECTSURF3D detects SURF points in a volume.
%
% INPUTS
% V:
%
% OUTPUTS
% points: a Nx5 array of detected points. Each point is represented as a
% row vector [x,y,z,scale,saliency].
%
% NOTES
% The SURF implementation in 2D expands the image 
% [citations here]
%% -----------------------------------------------------------------------
%V = 0;
koctaves = 2;
octaves = {[9,15,21,27],[15,27,39,51],[27,51,75,99],[51,99,147,195]};
J = integralimage3D(V);
[x0,y0,z0] = size(V);

% Calculate the responses from the Hessian based detector at all the scales
% and store the results in cell R.
R = cell(9,1);
for i = 1:koctaves
    fprintf(1, 'Calculating Hessians for scales...');
for scale = octaves{i}
    % Check to see if this scale has already been computed.
    if(length(R) < scale || isempty(R{scale}))
       fprintf(1, ' %i', scale);
       R{scale} = makedetH(J, scale);
    else
       fprintf(1, ' X');
    end
end
fprintf(' DONE.\n');
end

% Perform non-maxima supression using a 3x3x3x3 bounding box

% Interpolate responses in scale space in order to find response maxima.






end
