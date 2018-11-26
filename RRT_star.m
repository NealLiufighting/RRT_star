function RRT_star

    figure
    clf
    axes('box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
    hold on
    axis([0 1 0 1])
    daspect([1 1 1])
    title('RRT')    
       
    delay=0.02; % Rendering delay 

    C_init = point(0.2, 0.2); % Initial position
    C_goal = point(0.8, 0.8); % Target point
    N_steps=1000;
    p = 0.001;
    eps = 0.03; % Accuracy
    
    T = struct('vertex',[],'edge',[]);
     
    plot(C_init.x,C_init.y,'.','markersize',24,'color','b'); % init
    plot(C_goal.x,C_goal.y,'.','markersize',24,'color','r'); % goal
    obstacle=[point(0.3,0.3), point(0.8,0.3), point(0.7,0.7), point(0.3,0.7)];
    xs=[obstacle(1).x obstacle(2).x obstacle(3).x obstacle(4).x obstacle(1).x];
    ys=[obstacle(1).y obstacle(2).y obstacle(3).y obstacle(4).y obstacle(1).y];
    plot(xs,ys,'linewidth',2,'color','k'); % obstacle     
    
    pause(1);
    
    T.vertex=[C_init];    
    step=0;
    tree_size=1;
    parents = -1;
    costs = 0;    
    linehandles = plot(C_init.x,C_init.y);
    while(step<N_steps)
        
        C_rand=GenerateState();
        
        % Nearest tree vertex
        [C_near,index]=NearestNeighbour(C_rand,T); 
        
        h1 = plot([C_near.x C_rand.x], [C_near.y C_rand.y], 'k:');
        h2 = plot(C_rand.x, C_rand.y, '.k', 'markersize', 10);
        pause(delay);
        delete(h1);
        delete(h2);
        
        % Finding a conflict-free configuration on a segment
        C_new=FindStoppingState(C_near,C_rand,p,obstacle);  
        if(~isequaln(C_new,C_near))            
            
            tree_size=tree_size+1;           
            
            prel_handle = plot([C_near.x C_new.x],[C_near.y C_new.y],'b','linestyle','--');
            plot(C_new.x,C_new.y,'.k','markersize',10);              
            
            r=SearchRadius(step,N_steps); 
            
            hcircle = viscircles(gca,[C_new.x C_new.y],r,'color',[0.7 0.7 0.7],'linestyle','--');
                        
            % Search tree nodes inside the sphere of radius r with the center in C_new
            [C_Nearest_ind, distance_to_new]=NearestNeighbours(C_new,T,r); 
            if(~isempty(C_Nearest_ind))
            
                % Sorting of found vertexes at cost
                [sorted_costs,i_sort] = SortNearestNeighbours(C_Nearest_ind,costs,distance_to_new); 
                        
                % Search among the found vertexes available with minimum cost    
                [new_cost,i_parent] = MinCostParent(T,C_new,C_Nearest_ind,p,obstacle,sorted_costs,i_sort); 
                                        
                pause(delay);
                delete(prel_handle);          
                pause(delay);
                h_new = plot([T.vertex(C_Nearest_ind(i_parent)).x C_new.x],[T.vertex(C_Nearest_ind(i_parent)).y C_new.y],'b');
                pause(delay);
           
                T.vertex=[T.vertex, C_new];            
                parents = [parents C_Nearest_ind(i_parent)];
                costs = [costs new_cost];
            
                linehandles = [linehandles h_new];            
            
            % Optimization of paths inside the sphere with a new vertex
            [linehandles,parents,costs] = Rewire(T,tree_size,new_cost,distance_to_new,costs,parents,C_Nearest_ind,C_new,p,obstacle,linehandles,delay); 
            
            end
            
            delete(hcircle);
            
            success=(Distance(C_new,C_goal)<eps);
            
            if(success)
                OptimalPathPlot(T,parents);
                break;
            end
            
        end
        step=step+1;
    end
end

function res=point(x,y)
    res=struct('x',x,'y',y);
end

function res=line(p1,p2)
    res=struct('p1',p1,'p2',p2);
end

function res = GenerateState()
    x=rand;
    y=rand;
    res=point(x,y);
end

function [res,ind] = NearestNeighbour(P_rand, Tree)
    index=1;
    min = Distance(P_rand,Tree.vertex(1));
    for i=2:length(Tree.vertex)
        d=Distance(P_rand,Tree.vertex(i));
        if(min>d)
            min=d;
            index=i;
        end
    end 
    res=Tree.vertex(index);
    ind=index;
end

function res = Distance(p1,p2)
    res = sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2);
end

function res = CollisionFree(P1,P2,p,obstacle)
    res=0;
    P=FindStoppingState(P1,P2,p,obstacle);
    if(isequaln(P,P2))
        res=1;
    end
end

function res = FindStoppingState(P_near, P_rand, p, obstacle)
    vec=point(P_rand.x-P_near.x,P_rand.y-P_near.y);
    len_vec=sqrt(vec.x^2+vec.y^2);
    k=len_vec/p;
    lines=[];
    for i=1:length(obstacle)
        if(i<length(obstacle))
            lines=[lines,line(obstacle(i),obstacle(i+1))];
        else
            lines=[lines,line(obstacle(i),obstacle(1))];
        end
    end
    cross=0;
    res=P_near;
    for i=1:fix(k)
        P_test=point(P_near.x+i*(vec.x/k),P_near.y+i*(vec.y/k));
        L_test=line(P_near,P_test);
        for j=1:length(lines)
            cross=IsCross(L_test,lines(j));
            if(cross)
                break;
            end
        end
        if(cross)
           break;
        else
           res=P_test;
        end
    end
    if(~cross)
        res=P_rand;
    end
end

function res = IsCross(L1,L2)
    x=[L1.p1.x L1.p2.x L2.p1.x L2.p2.x]; 
    y=[L1.p1.y L1.p2.y L2.p1.y L2.p2.y];  
    dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]); 
    dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]); 

    if(dt1<=0 && dt2<=0) 
        res=1; 
    else 
        res=0; 
    end 
    
end

function res = SearchRadius(V,N_steps)
%     YRRG=2*(1+1/2)^(1/2)*(1/1);
%     res=YRRG*(log(V)/V)^(1/2);
    res=1*(1-V/N_steps);
end

function [inds,dist_to_new] = NearestNeighbours(P_new, Tree, radius)
    
    inds=[];
    dist_to_new=[];
    for i=1:length(Tree.vertex)
        dist=Distance(P_new,Tree.vertex(i));
        if(dist<radius)
            inds=[inds,i];
        end
        dist_to_new=[dist_to_new,dist];
    end
end

function [sorted_costs,i_sort] = SortNearestNeighbours(C_Nearest_ind,costs,distance_to_new) 
    [sorted_costs,i_sort] = sort(costs(C_Nearest_ind) + distance_to_new(C_Nearest_ind));
end

function [new_cost,i_parent] = MinCostParent(T,C_new,C_Nearest_ind,p,obstacle,sorted_costs,i_sort)
    found_collision_free_parent = false;
    i = 1;
    while ~found_collision_free_parent
        i_parent = i_sort(i);
        if (CollisionFree(T.vertex(C_Nearest_ind(i_parent)),C_new,p,obstacle))
            new_cost = sorted_costs(i);
            found_collision_free_parent = true;
        else
            i = i + 1;
        end
    end
end

function [linehandles,parents,costs] = Rewire(T,tree_size,new_cost,distance_to_new,costs,parents,C_Nearest_ind,C_new,p,obstacle,linehandles,delay)
    neighbors_to_rewire = find(new_cost + distance_to_new(C_Nearest_ind) < costs(C_Nearest_ind));
    for ir = 1:length(neighbors_to_rewire)
        i_node = C_Nearest_ind(neighbors_to_rewire(ir));
        if (CollisionFree(T.vertex(i_node),C_new,p,obstacle))
            parents(i_node) = tree_size;
            costs(i_node) = new_cost + distance_to_new(i_node);
            
            prel_handle =  plot([T.vertex(i_node).x C_new.x],[T.vertex(i_node).y C_new.y],'m');  
            pause(delay);
            delete(linehandles(i_node));
            pause(delay);
            delete(prel_handle);
            linehandles(i_node) = plot([T.vertex(i_node).x C_new.x],[T.vertex(i_node).y C_new.y],'b');
            pause(delay);
        end
    end
end

function OptimalPathPlot(T, parents)
    ix = length(T.vertex);
    while (parents(ix) ~= -1)
        plot([T.vertex(parents(ix)).x T.vertex(ix).x],[T.vertex(parents(ix)).y T.vertex(ix).y],'color',[0 0 0],'linewidth',2)
        plot([T.vertex(parents(ix)).x T.vertex(ix).x],[T.vertex(parents(ix)).y T.vertex(ix).y],'color',[0 0.9 0],'linewidth',1.5)
        ix = parents(ix);
    end
end