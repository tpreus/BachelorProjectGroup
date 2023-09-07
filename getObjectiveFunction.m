%function to calculate the objective vector c 
function objective = getObjectiveFunction(VX, FX, FeatX, VY, FY, FeatY)

    nx = size(VX,1);
    ny = size(VY,1);
    fx = size(FX,1);
    fy = size(FY,1);
    
    
    TR = triangulation(FY, VY);
    EY = edges(TR);
    my = size(EY,1);
    
    TR = triangulation(FX, VX);
    EX = edges(TR);
    mx = size(EX,1);

    objective = zeros(3*fx*ny + fx*ny, 1);
   
    for i = 1:ny
        for j = 1:fx
            objective( 3*(i-1)* fx + 3*(j-1)+1, 1) = norm(FeatX(FX(j,1), :)- FeatY(i, :) , 1);
            objective( 3*(i-1)* fx + 3*(j-1)+2, 1) = norm(FeatX(FX(j,2), :)- FeatY(i, :) , 1);
            objective( 3*(i-1)* fx + 3*(j-1)+3, 1) = norm(FeatX(FX(j,3), :)- FeatY(i, :) , 1);
        end
    end

    tmp = sum(objective) / (3*fx*ny);
    objective = sparse([objective; zeros(2*mx*my,1) + tmp]);% zeros(2*my*fx,1)+ sum(objective)]);

    disp("got objective function")

end