function ret = getCoordinates(VX, FX, VY, FY, x)

    function [pointX, pointY, pointZ] = getPoint(d1, d2, d3, v1, v2 ,v3) %return 3D coordinates
        disp("got point with weight " + d1 + " " + d2 + " " + d3 );
        pointX = d1 * VX(v1,1) + d2* VX(v2,1) + d3*VX(v3,1);
        pointY = d1 * VX(v1,2) + d2* VX(v2,2) + d3*VX(v3,2);
        pointZ = d1 * VX(v1,3) + d2* VX(v2,3) + d3*VX(v3,3);
    end




    nx = size(VX,1);
    fx = size(FX,1);
    ny = size(VY,1);
    xL = x(1:3*fx*ny,1);
    xI = x(3*fx*ny+1 : 3*fx*ny+fx*ny ,1);

    ret = zeros(ny, 3);

    for i = 1:ny %point looking at
        for j = 1:fx % triangle looking at
            if xI( (i-1)*fx + j) == 1 % point i in the j-th triangle 
                tmp = (3*fx*(i-1) + 3*(j-1) +1);

                d1 = xL(tmp);
                d2 = xL(tmp+1);
                d3 = xL(tmp+2);
    
                [p1, p2, p3] = getPoint(d1, d2 , d3, FX(j,1), FX(j,2), FX(j,3) );
                ret(i, 1) = p1;
                ret(i, 2) = p2;
                ret(i, 3) = p3;
            end
        end
    end

end