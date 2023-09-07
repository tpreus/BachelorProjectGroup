
forceEqualTriangulationOnDownsampledShapes = true;
experimentName = "test";
Xexperiment = ["100Ftoy_1.ply", "200Fsphere.ply", "200Fsphere.ply", "200Fsphere.ply", "200Fsphere.ply", "200Fsphere.ply" ,"200Ftoy_1.ply"];
Yexperiment = ["100Ftoy_2.ply", "200Fapple.ply", "200Fapple.ply", "200Fapple.ply", "200Fapple.ply", "200Fapple.ply", "toy_2.ply"];
numberFeatures = 50;



num_facesX = 60;
num_facesY = 60;
runtimes = zeros(size(numberOfTri,2),2);


% A_leq3, A_leq4, A_leq5, A_leq6
ListConstrains =[3,1,1,0];
    
solver = "gurobi"; % solver in {"gurobi", "mosek"}; %setting solver
timeout = 0.8; % hour
relGap = 1e-1;
nThreads = 4;

showResults = true;


for i = 1:1%size(numberOfTri,2)
    try 
        constrains = ListConstrains(i,:);
        num_facesX = numberOfTri(1);
        num_facesY = numberOfTri(1);
        [VX, FX] = readPLY(Xexperiment(i));
        [VY, FY] = readPLY(Yexperiment(i));
        
    
        [VX, FX, ~, IX] = decimate_libigl(VX, FX, num_facesX, 'Method', 'qslim');
        [VY, FY, ~, IY] = decimate_libigl(VY, FY, num_facesY, 'Method', 'qslim');
        
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
    
        FeatX = calcWks(VX, FX, numberFeatures, 6, 25);
        FeatY = calcWks(VY, FY, numberFeatures, 6, 25);
    
        nx = size(VX,1);
        ny = size(VY,1);
        fx = size(FX,1);
        fy = size(FY,1);
    
        [A_equal, A_leq, b_equal, b_leq] = getConstrains(VX, FX, VY, FY, constrains);
 
        
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
        
        xL = sdpvar(3*fx*ny, 1);
        xI = binvar(fx*ny,1); %already in \{0,1\}
        
        gamma = sdpvar(2*my*mx,1);
        constraint_values_in_zero_one = [0 <= xL; xL <= 1; gamma >= 0];
        % 
        x = [xL; xI; gamma];
        c = getObjectiveFunction(VX, FX, FeatX, VY, FY, FeatY);
        if ListConstrains(i,1) ==3
            sigma = sdpvar(2*my*fx,1);
            x = [x; sigma];
            c = [c; zeros(2*my*fx,1)+ sum(c)];
            constraint_values_in_zero_one = [constraint_values_in_zero_one; sigma >= 0];
        
        end

        objective = c'*x;
        
        disp(size(A_leq));
        disp(size(A_equal));
        disp(size(x));
        disp(size(b_leq));
        disp(size(b_equal));
        constraints = [ A_leq*x <= b_leq; A_equal*x == b_equal;  constraint_values_in_zero_one];
        % 
        
        
        disp("start optimization");
        tic;
        yalmipOut = optimize(constraints, objective, opts);
        t = toc;
        runtimes(i,1) = t;
        runtimes(i,2) = numberOfTri(i);
        XValue = value(x);
        
        points = getCoordinates(VX, FX, VY, FY, XValue);
        
        xL = XValue(1:3*fx*ny, 1);
        xI = XValue(3*fx*ny+1 : 3*fx*ny+fx*ny ,1);
        xGamma = XValue(3*fx*ny+fx*ny: 3*fx*ny+fx*ny+2*mx*my ,1);
        
        

        if ListConstrains(i,1) == 3
            xSigma = XValue(3*fx*ny+fx*ny+2*mx*my +1 : 3*fx*ny+fx*ny+2*mx*my + 2*fx*ny);
            cSigma = c(3*fx*ny+fx*ny+2*mx*my +1 : 3*fx*ny+fx*ny+2*mx*my + 2*fx*ny);
        end
        cL = c(1:3*fx*ny, 1);
        cI = c(3*fx*ny+1 : 3*fx*ny+fx*ny ,1);
        cGamma = c(3*fx*ny+fx*ny: 3*fx*ny+fx*ny+2*mx*my ,1);
        color = zeros(size(VY,1),3);
        
        showResults = true;
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
            savefig((string(i)+experimentName+ 'Fig.fig'));
        end
        disp("experiment " + string(i) + " finished successfully");
    catch
      disp("experiment " + string(i) + " failed");
    end

end


