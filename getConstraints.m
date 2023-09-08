function [A_equal, A_leq, b_equal, b_leq] = getConstraints(VX, FX, VY, FY, constraints)
    function ret = euclideanDist(i, j, V)
        ret = sqrt( (V(i,1) - V(j,1))^2 + (V(i,2) - V(j,2))^2 + (V(i,3) - V(j,3))^2 );
    end
    disp("start getting constraints");
    
%x will be (d_1^{v_1}, d_2^{v_1}, ... d_{n_x}^{v_1} , d_1^{v_2}, ..., d_{n_x}^{v_2}, ... , d_{n_y}^{v_{1}} ,..., d_{n_y}^{v_{v_{n_x} , )
    
    nx = size(VX,1);
    ny = size(VY,1);
    fx = size(FX,1);
    fy = size(FY,1);
    
    
    %--------
    %define A_eqal (The essential constraints)
    %--------
    % first \sum_{i=1}^3 d_i^{pf} = P_f^p
    AL1 = sparse(fx * ny, 3*fx*ny);
    b1 = sparse(fx*ny,1);

    %A2 contains \sum_{f \in F_x} P_p^f = 1
    A2 = sparse(ny, 3*fx*ny + fx*ny);
    b2 = ones(ny,1);
    for i = 1:ny*fx
        AL1(i, 3*(i-1)+1 : 3*(i-1) +3 ) = 1;
        if i <= ny
            A2(i, 3*fx*ny + (i-1)*fx +1 : 3*fx*ny + i*fx) = 1;
        end
    end
    
    AI1 = -speye(fx*ny);
    A1 = [AL1, AI1];
    A1 = sparse(A1);

    
    %--------
    %define A_leq (all optional constraints)
    %-------- 
    
    %Condition that all neighboring points from Y also must be in
    %neighboring triangles on X
    TR = triangulation(FY, VY);
    EY = edges(TR);
    my = size(EY,1);

    TR = triangulation(FX, VX);
    EX = edges(TR);
    mx = size(EX,1);

    A_leq3_1 = sparse(2*my*fx ,fx*ny); %constraints for directeley neighboring triangles
    b_leq3_1 = sparse(2*my*fx,1);

    if constraints(1) == 1
        for i = 1:2*my
            if i <= my
                v = EY(i, 1);
                w = EY(i, 2);
            elseif i > my
                w = EY(i-my, 1);
                v = EY(i-my, 2);
            end
    
            for j = 1:fx 
                A_leq3_1( (i-1) * fx + j, (v-1)*fx+ j) = 1;%if v is placed in triangle j
                c1 = FX(j,1);
                c2 = FX(j,2);
                c3 = FX(j,3);
                for k = 1:fx %find all neigboring vertices
                    %condition that vertex is in triangle that neighbours the
                    %triangle of v indirectly
                    if (c1 == FX(k,1) || c1 == FX(k,2) || c1 == FX(k,3) || c2 == FX(k,1) || c2 == FX(k,2) || c2 == FX(k,3) || c3 == FX(k,1) || c3 == FX(k,2) || c3 == FX(k,3) )
    %                 if (c1 == k1 && (c2 == k2 || c2 == k3)) || (c1 == k2 && (c2 == k3 || c2 == k1)) || (c1 == k3 && (c2 == k1 || c2 == k2)) || (c2 == k1 && (c3 == k2 || c3 == k3)) || (c2 == k2 && (c3 == k1 || c3 == k3)) || (c2 == k3 && (c3 == k1 || c3 == k2)) || (c3 == k1 && (c1 == k2 || c1 == k3)) || (c3 == k2 && (c1 == k1 || c1 == k3)) || (c3 == k3 && (c1 == k1 || c1 == k2))
                        A_leq3_1((i-1) * fx + j, (w-1)*fx + k) = -1;
                    end
                end
            end
    
        end

        A_leq3_1 = sparse([sparse(2*my*fx, 3*fx*ny) , A_leq3_1, sparse(2*my*fx, 2*mx*my)]);%, -speye(2*my*fx, 2*my*fx)]);
    
        A_leq3 = A_leq3_1;
        b_leq3 = b_leq3_1;
        disp("finished A_leq3");
    end



    if constraints(1) ==2 || constraints(1) ==3
        A_leq3_2 = sparse(2*my*fx ,fx*ny); %constraints for directeley neighboring triangles
        b_leq3_2 = sparse(2*my*fx,1);
        A_leq3_1 = sparse(2*my*fx ,fx*ny); %constrains for directeley neighboring triangles
        b_leq3_1 = sparse(2*my*fx,1);

        for i = 1:2*my
            if i <= my
                v = EY(i, 1);
                w = EY(i, 2);
            elseif i > my
                w = EY(i-my, 1);
                v = EY(i-my, 2);
            end
    
            for j = 1:fx 
                A_leq3_1( (i-1) * fx + j, (v-1)*fx+ j) = 1;%if v is placed in triangle j 
                A_leq3_2( (i-1) * fx + j, (v-1)*fx+ j) = 1;
    
                c1 = FX(j,1);
                c2 = FX(j,2);
                c3 = FX(j,3);
                for k = 1:fx %find all neigboring vertices
                    %condition that vertex is in triangle that neighbours the
                    %triangle of v indirectly
                    k1 = FX(k,1);
                    k2 = FX(k,2);
                    k3 = FX(k,3);
                    if (c1 == FX(k,1) || c1 == FX(k,2) || c1 == FX(k,3) || c2 == FX(k,1) || c2 == FX(k,2) || c2 == FX(k,3) || c3 == FX(k,1) || c3 == FX(k,2) || c3 == FX(k,3) )
                       A_leq3_1((i-1) * fx + j, (w-1)*fx + k) = -1;
                    %condition that vertex is in triangle that neighbours v
                    %directly
                        A_leq3_2((i-1) * fx + j, (w-1)*fx + k) = -1; 
                        for l = 1:fx
                            l1 = FX(l,1);
                            l2 = FX(l,2);
                            l3 = FX(l,3);
                            if(k1 == FX(l,1) ||k1 == FX(l,2) || k1 == FX(l,3) || k2 == FX(l,1) || k2 == FX(l,2) || k2 == FX(l,3) || k3 == FX(l,1) || k3 == FX(l,2) || k3 == FX(l,3) )
                               A_leq3_2((i-1) * fx + j, (w-1)*fx + l) = -1;
                            end
                        end
                    end
                end
            end
        end

        A_leq3_1 = sparse([sparse(2*my*fx, 3*fx*ny) , A_leq3_1, sparse(2*my*fx, 2*mx*my)]);%, -speye(2*my*fx, 2*my*fx)]);
        A_leq3_2 = sparse([sparse(2*my*fx, 3*fx*ny) , A_leq3_2, sparse(2*my*fx, 2*mx*my)]);%, sparse(2*my*fx, 2*my*fx)]);
    
        disp("finished A_leq3 2");
    end

    %------



    %condition that in every patch of triangles must be at least one vertex
    if constraints(3) == 1
        A_leq5 = sparse(fx ,fx*ny);
        for i = 1:ny
            for j = 1:fx
                c1 = FX(j,1);
                c2 = FX(j,2);
                c3 = FX(j,3);
                for k = 1:fx
                    k1 = FX(k,1);
                    k2 = FX(k,2);
                    k3 = FX(k,3);
                    if (c1 == FX(k,1) || c1 == FX(k,2) || c1 == FX(k,3) || c2 == FX(k,1) || c2 == FX(k,2) || c2 == FX(k,3) || c3 == FX(k,1) || c3 == FX(k,2) || c3 == FX(k,3) )
                        A_leq5(j, (i-1)*fx + k) = -1;
                    end
                end
            end
        end
        b_leq5 = sparse(fx,1) -1;

        A_leq5 = [sparse(fx, 3*fx*ny), A_leq5 , sparse(fx, 2*mx*my)];
        disp("finished A_leq5");
    end
    %-----

    %constraint for the triangle inequality
    if constraints(2) == 1
        A_leq4 = sparse(2*my * mx , 3*fx *ny + fx*ny);
        b_leq4 = zeros(2*my * mx , 1); 
        for e = 1:my% betrachte alle Kanten von Y
            v = EY(e,1);
            w = EY(e,2); 
            e_dist = euclideanDist(v,w, VY);
            for i = 1:2*mx
                if i <= mx 
                    c1 = EX(i,1);
                    c2 = EX(i,2);
                elseif i <= 2*mx
                    c1 = EX(i -mx,2);
                    c2 = EX(i-mx,1);
                elseif i <= 2*mx + vx
                    c1 = VX(i-2*mx);
                    c2 = VX(i-2*mx);
                end

                %find the two triangles that contain both c1, c2
                %\{c1, c2\} ist die Kante die die beiden Dreiecke, die gerade betrachtet werden trennt
                
    
                for k = 1:fx %finde die genauen Dreiecke, in die v und w gelegt werden soll und die von {c1, c2} getrennt werden
                    if FX(k,1) == c1 && (FX(k,2) == c2 || FX(k,3) == c2) || (FX(k,2) == c1 && (FX(k,1) == c2 || FX(k,3) == c2)) || (FX(k,3) == c1 && (FX(k,1) == c2 || FX(k,2) == c2))
                        if FX(k,1) == c1 %find the node of the triangle corresponding to c1 
                            vertexInTriangle = 1;
                        elseif(FX(k,2) == c1)
                            vertexInTriangle = 2;
                        elseif(FX(k,3) == c1)
                            vertexInTriangle = 3;
                        end
    
                        %approximate e_dist a by assuming that both triangles
                        %are equilateral
                        %approximate the edge length by the average of the edge
                        %length of the triangle
                        averageEdgeLength = ( euclideanDist( FX(k,1), FX(k,2), VX) + euclideanDist(FX(k,2), FX(k,3), VX) + euclideanDist(FX(k,3), FX(k,1), VX) )/3;
                        e_dist1 = (averageEdgeLength^2 + (averageEdgeLength/2)^2 )^0.5;
    
                        %both edges of the triangle k that are adjacent to c1
                        vertexInTriangle = vertexInTriangle -1; %wird kurzzeitig -1 gesetz, da dann mod besser berechenbar 
                        edgeC1W = euclideanDist( FX(k, vertexInTriangle +1), FX(k, mod(vertexInTriangle +1,3)+1 ), VX );
                        edgeC1U = euclideanDist( FX(k, vertexInTriangle +1), FX(k, mod(vertexInTriangle +2, 3)+1 ), VX );
                        averageLength = (edgeC1W + edgeC1U)/2;
                        vertexInTriangle = vertexInTriangle +1;
    
    
                        A_leq4((e-1)*2*mx + i , (v-1)*3*fx + 3*(k-1) +  vertexInTriangle) = averageLength; %d_i^{v f_k} = 1 %weight of vertex v in triangle k
                        b_leq4((e-1)*2*mx +i) =+ averageLength;
    
                        A_leq4((e-1)*2*mx +i, (3*fx+ny) + (v-1)*fx + k) = e_dist; 
                    end
                end
            end
        end
        disp("finished A_leq4");
        A_leq4 = [A_leq4, -speye(2*my * mx, 2*mx*my)];
    end

%---------------

%condition that at most a fixed number of points can lie in one triangle
    if constraints(4) == 1
        A_leq6 = sparse(fx, fx*ny);
        b_leq6 = sparse(fx ,1) + 1;
        for i = 1:fx
            for j = 1:ny
                A_leq6(i, (j-1)*fx +i ) = 1;
            end
        end
        A_leq6 = sparse([sparse(fx, 3*fx*ny) , A_leq6 , sparse(fx,2*mx*my)]);%, sparse(fx, 2*my*fx);]);
        disp("finished A_leq6");
    end
%-----

%merging the different copnstraints
    A_leq = sparse(1, 3*fx*ny + fx*ny + 2*mx*my);
    b_leq = sparse(1,1);
    A_equal = sparse([A1 ; A2]);
    A_equal = [A_equal, sparse(size(A_equal,1), 2*my*mx )];%2*my*fx
    b_equal = [b1 ; b2];
    if constraints(1) == 3
        A_equal = [A_equal, sparse(size(A_equal,1), 2*my*fx)];
        A_leq = [A_leq, sparse(1, 2*my*fx)];
        A_leq = [A_leq; [A_leq3_1, -speye(2*my*fx, 2*my*fx)]; [A_leq3_2, sparse(2*my*fx,2*my*fx) ]];
        b_leq = [b_leq; b_leq3_1; b_leq3_2];


        if constraints(2) == 1
            A_leq = [A_leq; [A_leq4, sparse(2*my*mx, 2*my*fx)]];
            b_leq = [b_leq; b_leq4];
        end
        if constraints(3) == 1
            A_leq = [A_leq; [A_leq5,sparse(fx, 2*my*fx)]];
            b_leq = [b_leq; b_leq5];
        end
        if constraints(4) == 1
            A_leq = [A_leq; [A_leq6, sparse(fx, 2*my*fx)]];
            b_leq = [b_leq; b_leq6];
        end
    else
        if constraints(1) == 1
            A_leq = [A_leq; A_leq3_1];
            b_leq = [b_leq; b_leq3_1];
        end
        if constraints(1) == 2
            A_leq = [A_leq; A_leq3_2];
            b_leq = [b_leq; b_leq3_2];
        end
        if constraints(2) == 1
            A_leq = [A_leq; A_leq4];
            b_leq = [b_leq; b_leq4];
        end
        if constraints(3) == 1
            A_leq = [A_leq; A_leq5];
            b_leq = [b_leq; b_leq5];
        end
        if constraints(4) == 1
            disp(size(A_leq));
            disp(size(A_leq6));
            A_leq = [A_leq; A_leq6];
            b_leq = [b_leq; b_leq6];
        end
    end

    
    
    disp("got all constraints");
end