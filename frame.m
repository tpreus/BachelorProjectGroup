%frame for the optimization problem, allows to compute multiple instances


forceEqualTriangulationOnDownsampledShapes = true;
experimentName = "test";
Xexperiment = ["100Fsphere.ply"];
Yexperiment = ["100Fapple.ply"];
numberFeatures = 50;

num_facesX = 200;
num_facesY = 200;

runtimes = zeros(size(numberOfTri,2),2);


%setting the constraints that will be used:
%first entry is the geometric consistency constraint with 0 for not using
%it, 1 for the normal constraint, 2 for the double neighborhood and 2 for
%the combination of single and double neighborhood with relaxation
%variables
%second entry is for the triangle inequality \in \{0,1\}
%third is for the lower bound on numbers of points per triangle patch \in
%\{0,1\}
%fourth is for the upper bound on the number of points per single triangle
ListConstraints =[3,1,1,0];
    
solver = "gurobi"; % solver in {"gurobi", "mosek"}; %setting solver
timeout = 0.8; % hour
relGap = 1e-1;
nThreads = 4;

showResults = true;


for i = 1:1%
    % try 
        constraints = ListConstraints(i,:);

        [VX, FX] = readPLY(Xexperiment(i));
        [VY, FY] = readPLY(Yexperiment(i));
        
    
        % [VX, FX, ~, IX] = decimate_libigl(VX, FX, num_facesX, 'Method', 'qslim');
        % [VY, FY, ~, IY] = decimate_libigl(VY, FY, num_facesY, 'Method', 'qslim');
        % 
        %Verschieben ------
        tmp = min(min(VY));
        VY = VY - tmp;  
        tmp = max(max(VY));
        VY = VY/tmp;
        % 
        tmp = min(min(VX));
        VX = VX - tmp;  
        tmp = max(max(VX));
        VX = VX/tmp;
        VX = VX - mean(VX);
        VY = VY - mean(VY);
        %------
    
        FeatX = calcWks(VX, FX, numberFeatures, 6, 25);%calculating the feature vectors
        FeatY = calcWks(VY, FY, numberFeatures, 6, 25);
    
        nx = size(VX,1);
        ny = size(VY,1);
        fx = size(FX,1);
        fy = size(FY,1);
    
        [A_equal, A_leq, b_equal, b_leq] = getConstraints(VX, FX, VY, FY, constraints); %calculating constraints for the problem
 
        
        if solver == "mosek"
            opts = sdpsettings('solver', 'mosek','verbose', 1);
            opts.mosek.MSK_DPAR_MIO_MAX_TIME = timeout*3600;
            opts.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = timeout*3600;
            opts.mosek.MSK_DPAR_MIO_TOL_REL_GAP = relGap;
            opts.mosek.MSK_IPAR_NUM_THREADS = nThreads;
        elseif  solver == "gurobi"
            opts = sdpsettings('solver','gurobi','verbose',1);
            opts.gurobi.TimeLimit = timeout*3600;
            opts.gurobi.MIPGap = relGap;
            opts.gurobi.Threads =nThreads;
        end
        
        TR = triangulation(FY, VY);
        EY = edges(TR);
        my = size(EY,1);
        
        TR = triangulation(FX, VX);
        EX = edges(TR);
        mx = size(EX,1);
        
        xL = sdpvar(3*fx*ny, 1); %linear optimization variables
        xI = binvar(fx*ny,1); % binary optimization variables already in \{0,1\}
        
        gamma = sdpvar(2*my*mx,1); %slack variable for the triangle inequality
        constraint_values_in_zero_one = [0 <= xL; xL <= 1; gamma >= 0];
        % 
        x = [xL; xI; gamma];
        c = getObjectiveFunction(VX, FX, FeatX, VY, FY, FeatY);
        if ListConstraints(i,1) ==3 % only relevant if constraint for double neighborhood with relaxation is added
            sigma = sdpvar(2*my*fx,1);
            x = [x; sigma];
            c = [c; zeros(2*my*fx,1)+ sum(c)];
            constraint_values_in_zero_one = [constraint_values_in_zero_one; sigma >= 0];
        
        end

        objective = c'*x;

        constraints = [ A_leq*x <= b_leq; A_equal*x == b_equal;  constraint_values_in_zero_one];

        disp("start optimization");
        tic;%to track runtime
        yalmipOut = optimize(constraints, objective, opts);
        t = toc;
        runtimes(i,1) = t;
        runtimes(i,2) = numberOfTri(i);
        
        
        XValue = value(x);
        xL = XValue(1:3*fx*ny, 1);
        xI = XValue(3*fx*ny+1 : 3*fx*ny+fx*ny ,1);
        xGamma = XValue(3*fx*ny+fx*ny: 3*fx*ny+fx*ny+2*mx*my ,1);

        points = getCoordinates(VX, FX, VY, FY, XValue); %calculating the points on mesh $X$

        if ListConstraints(i,1) == 3
            xSigma = XValue(3*fx*ny+fx*ny+2*mx*my +1 : 3*fx*ny+fx*ny+2*mx*my + 2*fx*ny);
            cSigma = c(3*fx*ny+fx*ny+2*mx*my +1 : 3*fx*ny+fx*ny+2*mx*my + 2*fx*ny);
        end
        cL = c(1:3*fx*ny, 1);
        cI = c(3*fx*ny+1 : 3*fx*ny+fx*ny ,1);
        cGamma = c(3*fx*ny+fx*ny: 3*fx*ny+fx*ny+2*mx*my ,1);
        color = zeros(size(VY,1),3);
        
        if showResults == true
            color(:,1) = (VY(:, 1)+ 0.001 - min(VY(:,1)))/(0.01+ max(VY(:,1)) - min(VY(:,1)));
            color(:,2) = (VY(:, 2)+ 0.001 - min(VY(:,2)))/(0.01+ max(VY(:,2)) - min(VY(:,2)));
            color(:,3) = (VY(:, 3)+ 0.001 - min(VY(:,3)))/(0.01+ max(VY(:,3)) - min(VY(:,3)));

            figure('units','normalized','position',[0,0,1,1]);
            set(gcf,'color','w');
            axis equal;
            subplot(1,2,1);
            axis equal
            axis off
            title("figure X with points from Y")
            p2 = patch('Faces', FX, 'Vertices', VX, 'FaceColor', 'blue');
            p2.FaceVertexAlphaData = 0;
            p2.FaceAlpha = 'flat'; 
            hold on;

            p1 = patch('Faces',FY,'Vertices',points,'FaceVertexCData',color,'FaceColor','interp');
            p1.FaceVertexAlphaData = 0.5;
            p1.FaceAlpha = 'flat'; 
            hold on
            scatter3(points(:,1), points(:,2), points(:,3), 70, color, 'filled');
        
            subplot(1,2,2);
            title("figure Y")
            patch('Faces',FY,'Vertices',VY,'FaceVertexCData',color,'FaceColor','interp');
           
            hold on
            h = scatter3(VY(:,1), VY(:,2), VY(:,3), 70, color, 'filled');
            hold off
            axis equal
            axis off
        end
        disp("experiment " + string(i) + " finished successfully");
    % catch
    %   disp("experiment " + string(i) + " failed");
    % end

end


