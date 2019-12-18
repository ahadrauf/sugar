function [sigma1,sigma2] = sigma(alpha,w0,u,k,a)

sigma1 = (3/8 * alpha * a^2 / w0) + (k^2 /(4 * w0^2 * a^2) - u^2)^(1/2);
sigma2 = (3/8 * alpha * a^2 / w0) - (k^2 /(4 * w0^2 * a^2) - u^2)^(1/2);

