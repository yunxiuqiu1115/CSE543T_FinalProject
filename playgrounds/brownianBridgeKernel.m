% brownian bridge kernel

function kernel = brownianBridgeKernel(x, y)
    kernel = min(x,y) -x*y;
end
