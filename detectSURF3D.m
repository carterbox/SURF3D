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

koctaves = 2; nbhood = [1,1,1,1];
octaves = {[15,21,27,33,39,45],[39,51,63,75,87,99],[27,51,75,99],[51,99,147,195]};
J = integralimage3D(V);

%% Calculate the responses from the Hessian based detector
R = cell(9,1); % Store the results in cell R.
for i = 1:koctaves
    fprintf(1, '\nCalculating Hessians for scales...');
for scale = octaves{i}
    % Check to see if this scale has already been computed.
    if(length(R) < scale || isempty(R{scale}))
       fprintf(1, ' %i', scale);
       temp = makedetH(J, scale);
       %load(sprintf('Rscale%02i.mat',scale),'temp');
       save(sprintf('Rscale%02i.mat',scale),'temp');
       R{scale} = temp;
    else
       fprintf(1, ' X');
    end
end
fprintf(' DONE.');
end

clear('J');
%% Perform non-maxima supression using a 2*nbhood+1 bounding box
fprintf(1, '\nSupressing non-maxima in octaves...');
maxima = cell(koctaves,1);
for i = 1:koctaves
    fprintf(1, '\n%i', i);
    % Concat all the detector responses from the octave
    A = R{octaves{i}(1)};
    for j = octaves{i}(2:end), A = cat(4,A,R{j}); end
    maxima{i} = nonmaximumsupression(A,nbhood);
    
    % Interpolate responses in scale space in order to find response maxima.
    fprintf(1, ' refining peaks... ');
    maxima{i} = refinepeaks(maxima{i},A,4,octaves{i});
    fprintf(' DONE.');
end

%% Combine all the results and put strongest responses first.
points = maxima{1};
for i = 2:length(maxima), points = cat(1,points,maxima{i}); end
points = sortrows(points, -5); % reponse magnitude is col 5.

fprintf(1,'\nSUCCCESS\n');
end
