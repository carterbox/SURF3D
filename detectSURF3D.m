function [points,scale,saliency,signofLaplacian] = detectSURF3D(V)
% DETECTSURF3D detects SURF points in a volume. #parallelchild
%
% INPUTS
% V: the input volume.
%
% OUTPUTS
% points: a Nx5 array of detected points. Each point is represented as a
% row vector [x,y,z,scale,saliency].
%
% NOTES
% This work is an extension of the work by Bay, Herbert, et al. "Speeded-up
% robust features (SURF)." Computer vision and image understanding 110.3 
% (2008): 346-359.
%
%% -----------------------------------------------------------------------
[x0,y0,z0] = size(V);
log = 1; % 1 is the default terminal.
koctaves = uint8(1); nbhood = [3;3;3;3];
octaves = {[9,15,21,27]};
scales = {[1.2,1.6,2.1,2.3]};
%[27,33,39,45],[39,51,63,75],[27,51,75,99],[51,99,147,195]};
threshold = 10000;

%% Calculate the responses from the Hessian based detector
R = cell(27,1); S = cell(27,1); % Store the results in cell R.
for i = 1:koctaves
    % Pad the array for the largest filter in the octave.
    buffer = repmat((max(octaves{i})-1)/2,[1,3]);
    padded_V = padarray(V,buffer,'replicate','both');
    
    J = integralimage3D(padded_V);
    
    fprintf(log, '\nCalculating Hessians for filters...');
    for filter_size = octaves{i}
        % Check to see if this filter size has already been computed.
        if(length(R) < filter_size || isempty(R{filter_size}))
            filename = sprintf('./Rfilter%02i.mat',filter_size);
            if(exist(filename, 'file') > 0)
                % Load the saved file.
                fprintf(log, ' L');
                load(filename,'detH', 'signL');
                R{filter_size} = detH;
                S{filter_size} = signL;
            else
                % Compute it from scratch.
                fprintf(log, ' %i', filter_size);
                [detH, signL] = makedetH(J, filter_size);
                % Unpad the result.
                detH = detH(1+buffer:end-buffer,...
                            1+buffer:end-buffer,...
                            1+buffer:end-buffer);
                signL = signL(1+buffer:end-buffer,...
                    1+buffer:end-buffer,...
                    1+buffer:end-buffer);
                % Save it for later.
                save(filename,'detH','signL');
                R{filter_size} = detH;
                S{filter_size} = signL;
            end
        else
            % Already present.
            fprintf(log, ' X');
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

    % Concat all the detector responses from the octave.
    s0 = length(octaves{i});
    A(x0,y0,z0,s0) = single(0); 
    for j = 1:s0
        A(:,:,:,j) = R{octaves{i}(j)};
    end
    
    maxima{i} = nonmaximumsupression(A,nbhood);
    
    % Get sign values for each of the maxima
    B(x0,y0,z0,s0) = int8(0); 
    for j = 1:s0
        B(:,:,:,j) = S{octaves{i}(j)};
    end
    maxima{i} = getsigns(maxima{i},B);
    
    % Interpolate responses in scale space in order to find response maxima.
    fprintf(log, ' refining peaks... ');
    maxima{i} = refinepeaks(maxima{i},A,4,scales{i});
    
    fprintf(' DONE.');
end

% Combine all the results and put strongest responses first.
points = maxima{1};
for i = 2:length(maxima), points = cat(1,points,maxima{i}); end
points = sortrows(points, -5); % reponse magnitude is col 5.
numabovethresh = sum(points(:,5)>threshold);
points = points(1:numabovethresh,:);

scale = points(:,4);
saliency = points(:,5);
signofLaplacian = int8(points(:,6));
points = points(:,1:3);

fprintf(log,'\nSUCCCESS\n');

end

function [points] = getsigns(points,signgrid)

numpoints = size(points,1);
for i = 1:numpoints
    
    points(i,6) = signgrid(points(i,1),points(i,2),points(i,3),points(i,4));
    
end
end
