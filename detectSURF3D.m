function [points] = detectSURF3D(V)
% DETECTSURF3D detects SURF points in a volume. #parallelchild
%
% INPUTS
% V: the input volume.
%
% OUTPUTS
% points: a Nx5 array of detected points. Each point is represented as a
% row vector [x,y,z,e,saliency].
%
% NOTES
% This work is an extension of the work by Bay, Herbert, et al. "Speeded-up
% robust features (SURF)." Computer vision and image understanding 110.3 
% (2008): 346-359.
%
%% -----------------------------------------------------------------------
log = 1; % 1 is the default terminal.
koctaves = 1; nbhood = [3,3,3,3];
octaves = {[9,15,21,27]};
%[27,33,39,45],[39,51,63,75],[27,51,75,99],[51,99,147,195]};

%% Calculate the responses from the Hessian based detector
R = cell(9,1); % Store the results in cell R.
for i = 1:koctaves
    % Pad the array for the largest filter in the octave.
    % For an octive size of 5, is it more efficient to pad once or pad for
    % each filter separately?
    buffer = repmat((max(octaves{i})-1)/2,[1,3]);
    padded_V = padarray(V,buffer,'replicate','both');
    J = integralimage3D(padded_V);
    fprintf(log, '\nCalculating Hessians for filters...');
for filter_size = octaves{i}
    % Check to see if this filter size has already been computed.
    filename = sprintf('./Rfilter%02i.mat',filter_size);
    if(exist(filename, 'file') > 0)
        fprintf(log, ' L');
        load(filename,'detH');
        R{filter_size} = detH;
    else
        if(length(R) < filter_size || isempty(R{filter_size}))
           fprintf(log, ' %i', filter_size);
           padded_detH = makedetH(J, filter_size);
           % Unpad the result.
           detH = padded_detH(1+buffer:end-buffer,...
                              1+buffer:end-buffer,...
                              1+buffer:end-buffer);
           save(filename,'detH');
           R{filter_size} = detH;
        else
           fprintf(log, ' X');
        end
    end
end
fprintf(' DONE.');
end

clear('J');
%% Perform non-maxima supression using a 2*nbhood+1 bounding box
fprintf(log, '\nSupressing non-maxima in octaves...');
maxima = cell(koctaves,1);
for i = 1:koctaves
    fprintf(log, '\n%i nb:[%i,%i,%i,%i]', i, nbhood);
    % Concat all the detector responses from the octave
    A = R{octaves{i}(1)};
    for j = octaves{i}(2:end), A = cat(4,A,R{j}); end
    maxima{i} = nonmaximumsupression(A,nbhood);
    
    % Interpolate responses in scale space in order to find response maxima.
    fprintf(log, ' refining peaks... ');
    maxima{i} = refinepeaks(maxima{i},A,4,octaves{i});
    fprintf(' DONE.');
end
% Combine all the results and put strongest responses first.
points = maxima{1};
for i = 2:length(maxima), points = cat(1,points,maxima{i}); end
points = sortrows(points, -5); % reponse magnitude is col 5.

fprintf(log,'\nSUCCCESS\n');

end
