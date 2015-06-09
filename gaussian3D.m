function [g] = gaussian3D(partial) 
% GAUSSIAN3D creates a discretized and cropped 3D gaussian centered on a
% 9x9 grid.
% 
%% -----------------------------------------------------------------------
sigma = 1.2;
[X,Y,Z] = meshgrid(-4:4,-4:4,-4:4);

g = gaussianhelper(X,Y,Z, sigma, partial);

% Normalize the gaussian by the smallest and largest number
%big = max(max(max(g))); small = min(min(min(g)));
%g = (g - small)/(big-small);

end

function [t] = gaussianhelper(x,y,z, sigma, partial)

    switch partial
        case 0 % 3D gaussian
            t = 1/(sigma.^3 .* (2*pi).^(3/2)).* exp(-(x.^2 + y.^2 + z.^2)/sigma.^2);
        case 11 % Second partial derivative dxdx
            %  < t < 0.7792721
            t = (sqrt(2) .* x.^2 / (pi.^(3/2) .* sigma.^7) - 1 / (sqrt(2) .* pi.^(3/2) .* sigma.^5)).*exp(-(x.^2+y.^2+z.^2)/sigma.^2);
        case 12 % Second partial derivative dxdy
            % 0.50006 < t < 0.49994
            t = sqrt(2).*x.*y.*exp(-(x.^2+y.^2+z.^2)/sigma.^2)./(pi.^(3/2)*sigma.^7);
    end
    %t = 2.8284.*t;
end

