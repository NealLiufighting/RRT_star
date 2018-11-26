function RRT


    figure
    clf
    axes('box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
    hold on
    axis([0 1 0 1])
    daspect([1 1 1])
    title('RRT')    
    
    delay=0.5; % Задержка отрисовки

    C_init = point(0.2, 0.2); % Начальное положение
    C_goal = point(0.8, 0.8); % Целевая точка
    N_steps=1000;
    N_extend=10; % Попытка включить целевую точку в дерево
    p = 0.02; 
    
    T = struct('vertex',[],'edge',[]);
     
    plot(C_init.x,C_init.y,'.','markersize',24,'color','b'); % init
    plot(C_goal.x,C_goal.y,'.','markersize',24,'color','r'); % goal
    obstacle=[point(0.3,0.3), point(0.8,0.3), point(0.7,0.7), point(0.3,0.7)]; % препятствие
    xs=[obstacle(1).x obstacle(2).x obstacle(3).x obstacle(4).x obstacle(1).x];
    ys=[obstacle(1).y obstacle(2).y obstacle(3).y obstacle(4).y obstacle(1).y];
    plot(xs,ys,'linewidth',2,'color','k'); %      
    
    pause(1);
    
    T.vertex=[C_init];    
    step=0;
    success=0; % false
    tree_size=1;
    parents = -1;
    while((step<N_steps)&&(success==0))
        if(0 ~= mod(step,N_extend))
            C_rand=GenerateState();
        else
            C_rand=C_goal; % Попытка включить целевую точку в дерево
        end
        
        [C_near,index]=NearestNeighbour(C_rand,T); % Ближайшая вершина дерева
        
        h1 = plot([C_near.x C_rand.x], [C_near.y C_rand.y], 'k:');
        h2 = plot(C_rand.x, C_rand.y, '.k', 'markersize', 10);
        pause(delay);
        delete(h1);
        delete(h2);
        
        C_new=FindStoppingState(C_near,C_rand,p,obstacle); % Поиск безконфликтной конфигурации на отрезке
        
        if(~isequaln(C_new,C_near)) 
        
            plot([C_near.x C_new.x],[C_near.y C_new.y],'b');
            plot(C_new.x,C_new.y,'.k','markersize',10);
            pause(delay);
            
            T.vertex=[T.vertex, C_new];
            tree_size=tree_size+1;
            parents = [parents index];
            %T.edge=[T.edge, line(C_near,C_new)];
            success=(Distance(C_new,C_goal)<=p);
        end
        step=step+1;
    end
    if(success)
        OptimalPathPlot(T,parents);
    end
end

function res=point(x,y)
    res=struct('x',x,'y',y);
end

function res=line(p1,p2)
    res=struct('p1',p1,'p2',p2);
end

function res=edge(l,n1,n2)
    res=struct('line',l,'n1',n1,'n2',n2);
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

function OptimalPathPlot(T, parents)
    ix = length(T.vertex);
    while (parents(ix) ~= -1)
        plot([T.vertex(parents(ix)).x T.vertex(ix).x],[T.vertex(parents(ix)).y T.vertex(ix).y],'color',[0 0 0],'linewidth',2)
        plot([T.vertex(parents(ix)).x T.vertex(ix).x],[T.vertex(parents(ix)).y T.vertex(ix).y],'color',[0 0.9 0],'linewidth',1.5)
        ix = parents(ix);
    end
end