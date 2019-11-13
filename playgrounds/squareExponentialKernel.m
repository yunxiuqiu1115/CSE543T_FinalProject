% square exponential kernel

function kernel = squareExponentialKernel(x, y, lengthScale, magnitudeScale)
    kernel = magnitudeScale^2*exp(-(x-y)^2/2/lengthScale^2);
end
