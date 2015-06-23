function [g] = gaussian3D(fsize, varargin) 
% GAUSSIAN3D creates a discretized and cropped 3D gaussian centered on a
% square grid with side length FSIZE.
% gaussian3D(fsize)
% gaussian3D(fsize, sigma)
% gaussian3D(fsize, sigma, derivative)
%
% INPUTS
% fsize: scalar side length of the cropped gaussian filter. 
% sigma: the standard deviation of the gaussian.
% derivative: choosing 11 or 12 for this option returns either second
% partial derivative of the guassian instead of the non-differentiated
% gaussian.
%
% OUTPUTS
% g (single): a 3D matrix with the requested filter.
%
% NOTES
% The filter is normalized so the elements sum to 1.
% https://en.wikipedia.org/wiki/Multivariate_normal_distribution
%% -----------------------------------------------------------------------
% Set default values for optional inputs.
sigma0 = 1.2; derivative0 = 0; 

% Set input requirements.
p = inputParser;
addRequired(p,'fsize', @isnumeric);
addOptional(p,'sigma', sigma0, @isnumeric);
addOptional(p,'derivative', derivative0, @isnumeric);

% Assign values.
parse(p,fsize,varargin{:});
fsize  = single(p.Results.fsize);
sigma = single(p.Results.sigma);
derivative  = p.Results.derivative;
%% -----------------------------------------------------------------------

% Generate a list of points to evaluate the gaussian.
l0 = (fsize - 1)/2;
[X,Y,Z] = meshgrid(-l0:l0,-l0:l0,-l0:l0);

% Calculate values at each of the points.
g = gaussianhelper(X,Y,Z, sigma, derivative);

% Normalize the filter to a volume of 1
g = g/sum(g(:));

end

function [t] = gaussianhelper(x,y,z, sigma, derivative)
% GAUSSIANHELPER calculates the values of a 3D gaussian or its derivative
%   at point [x,y,z].
%% -----------------------------------------------------------------------
switch derivative
    case 0 % 3D gaussian
        t = exp(-(x.^2 + y.^2 + z.^2)/sigma.^2);
        %1/(sigma.^3 .* (2*pi).^(3/2)).*exp(-(x.^2 + y.^2 + z.^2)/sigma.^2);
    case 11 % Second partial derivative dxdx
        %t = (sqrt(2) .* x.^2 / (pi.^(3/2) .* sigma.^7) - 1 / (sqrt(2) .* pi.^(3/2) .* sigma.^5)).*exp(-(x.^2+y.^2+z.^2)/sigma.^2);
        t = (2.*(2*x.^2-sigma.^2)/sigma.^4).*exp(-(x.^2+y.^2+z.^2)/sigma.^2);
    case 12 % Second partial derivative dxdy
        %t = sqrt(2).*x.*y.*exp(-(x.^2+y.^2+z.^2)/sigma.^2)./(pi.^(3/2)*sigma.^7);
        t = (4.*x.*y./sigma.^4).*exp(-(x.^2+y.^2+z.^2)/sigma.^2);
end

end
