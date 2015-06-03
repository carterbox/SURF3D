function [g] = gaussian3D() 

sigma = 1.2;

[X,Y,Z] = meshgrid(-4:4,-4:4,-4:4);

g = gaussianhelperxx(X,Y,Z, sigma);
big = max(max(max(g))); small = min(min(min(g)));

g = (g - small)/(big-small);

end

function [t] = gaussianhelper(x,y,z, sigma)
	t = 1/(sigma.^3 .* (2*pi).^(3/2)).* exp(-(x.^2 + y.^2 + z.^2)/sigma.^2);
end

function [t] = gaussianhelperxx(x,y,z, sigma)

% t > 
% t < < 0.7792721
t = (sqrt(2) .* x.^2 / (pi.^(3/2) .* sigma.^7) - 1 / (sqrt(2) .* pi.^(3/2) .* sigma.^5)).*exp(-(x.^2+y.^2+z.^2)/sigma.^2);
end

function [t] = gaussianhelperxy(x,y,z,sigma)

% t > 0.50006
% t < 0.49994

t = sqrt(2).*x.*y.*exp(-(x.^2+y.^2+z.^2)/sigma.^2)./(pi.^(3/2)*sigma.^7);
end