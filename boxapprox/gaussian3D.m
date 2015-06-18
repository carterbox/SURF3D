function [g] = gaussian3D(size, varargin) 
% GAUSSIAN3D creates a discretized and cropped 3D gaussian centered on a
% square grid with side length SIZE.
% gaussian3D(size)
% gaussian3D(size, sigma)
% gaussian3D(size, derivative)
%% -----------------------------------------------------------------------
% Set default values for optional inputs.
sigma0 = 1.2; derivative0 = 0; 

% Set input requirements.
p = inputParser;
addRequired(p,'size', @isnumeric);
addOptional(p,'sigma', sigma0, @isnumeric);
addOptional(p,'derivative', derivative0, @isnumeric);

% Assign values.
parse(p,size,varargin{:});
size  = p.Results.size;
sigma = p.Results.sigma;
derivative  = p.Results.derivative;
%% -----------------------------------------------------------------------

% Generate a list of points to evaluate the gaussian.
l0 = (size - 1)/2;
[X,Y,Z] = meshgrid(-l0:l0,-l0:l0,-l0:l0);

% Calculate values at each of the points.
g = gaussianhelper(X,Y,Z, sigma, derivative);

% Normalize the gaussian by the smallest and largest number
%big = max(max(max(g))); small = min(min(min(g)));
%g = (g - small)/(big-small);

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
% t = 2.8284.*t; % Normalizes area to 1 for sigma = 1.2
end
