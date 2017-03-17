function intensity = GaussIntegrate(radius, w0, w1, w3, w5, w6)
c1 = -1/2*(1-w6)^2;
fun = @(x,y) w0 + w1.*exp(c1*((x./w3).^2+(y./w5).^2-(2*w6/w3/w5.*x.*y)));
intensity = integral2(fun,-radius,radius,-radius,radius);
end
