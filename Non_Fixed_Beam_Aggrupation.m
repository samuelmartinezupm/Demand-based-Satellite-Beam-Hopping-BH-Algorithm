function [S_interest]=Non_Fixed_Beam_Aggrupation(n_users, adj_matrix, algorithm_repetition)

S_interest=zeros(n_users); % In order to fulfill the condition in the first round the condition, S length is always goint to be lower or in the worst case equal to the user number.

for rep=1:algorithm_repetition   % To avoid local optima, the algorithm should be runned several times, and the best solution chosen (the one with in total lower number of cells).
    % [START] Algorithm
    
    % 0) Maximum cliques finding using Bron-Kerbosch algorithm:
    %A=[1 1 0 0 0 0; 1 1 1 1 0 0; 0 1 1 1 0 0; 0 1 1 1 0 1; 0 0 0 0 1 1; 0 0 0 1 1 1]
    %M=BK_MaxClique(A)
    maxcliques=BK_MaxClique(adj_matrix);
    
    % 1) Sort maximum cliques in descending order
    [val_maxcliques,idx_maxcliques]=sort_random_in_equality(sum(maxcliques));
    
    % 2) Inicialize emplty list that is used for saving the cliques that are chosen + define the collision (happens when two cliques share one or more vertex) degree (number of shared vertexes) variable:
    S={};
    collision=0;
    vertex_list=[];
    diff_user_count=0;
    count=1;
    
    % 3) Following the defined descending order, if the first selected clique
    % does not collide with any other clique in the solution set [S], it is
    % beign added into the the set [S]. Once the round has finished, the non collision criteria is increased till covering all the beams.
    
    while diff_user_count<n_users
        for i=1:size(maxcliques,2) % for each clique
            % The clique is: find(maxcliques(:,idx_maxcliques(i))) <- the non cero vertex indexes of the clique
            current_maxclique=find(maxcliques(:,idx_maxcliques(i)));
    
            if (length(intersect(current_maxclique,vertex_list))<= collision) % Add clique into the solution
                if ~isempty(setdiff(current_maxclique,vertex_list))
                    S{count}=setdiff(current_maxclique,vertex_list);
                    vertex_list=[vertex_list; S{count}];
                    diff_user_count=length(unique(vertex_list));
                    count=count+1;
                end
            end
    
        end
        collision=collision+1; % increase the collision allowance
    end
    % [END] Algorithm

    if length(S_interest)>=length(S)
        S_interest=S;
    end
end

end

function [val idx] = sort_random_in_equality(num_vertex)
[val,idx] = sort(num_vertex+(rand(size(num_vertex))-0.5),'descend'); % add a random value between +-0.5 so that in case of equality the random ordering is performed -> non deterministic algorithm
val = round(val); %integer return
end
   