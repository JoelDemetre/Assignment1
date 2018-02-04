function [ xinter1, xinter2, yinter1, yinter2 ] = BoxIntercept( x, y, xNext, yNext, xboxLow, xboxHigh, yboxLow, yboxHigh)  
xinter1 = 0;
xinter2 = 0;
yinter1 = 0;
yinter2 = 0;
m = (yNext - y)/(xNext-x);
b = y-m*x;

ycal = m*xboxLow + b;
if ycal < yboxHigh && ycal > yboxLow && ((ycal < y && ycal > yNext) || (ycal > y && ycal < yNext))
    xinter1 = 1;
    return;
end
ycal = m*xboxHigh + b;
if ycal < yboxHigh && ycal > yboxLow && ((ycal < y && ycal > yNext) || (ycal > y && ycal < yNext))
    xinter2 = 1;
    return;
end
xcal = (yboxLow-b)/m;
if xcal < xboxHigh && xcal > xboxLow && ((xcal < x && xcal > xNext) || (xcal > x && xcal < xNext))
    yinter1 = 1;
    return;
end
xcal = (yboxHigh-b)/m;
if xcal < xboxHigh && xcal > xboxLow && ((xcal < x && xcal > xNext) || (xcal > x && xcal < xNext))
    yinter2 = 1;
    return;
end

end

