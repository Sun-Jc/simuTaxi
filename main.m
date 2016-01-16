%% definitions:
% hypergraph:                                 nodes# - Nv
%                                             edges# - Ne
%                                   node_i_on_edge_j - G[Ne][Nv]
% conpensation:           node_i's_compensation      - C[Nv]
%                         total_compensation_on_e_j  - C_nom[Ne]
% possibility:         psb_that_node_i_served_by_e_j - P_ser[Ne][Nv]
%                      psb_that_node_i_get_served    - P_suc[Nv]
% utility:                             utility_for_i - U[Nv]
% const:    number_to_keep_denominator_larger_than_0 - alp
%                                             max_Nv - max_Nv
%                                             min_nv - min_Nv
%                                          max_Nv/Ne - max_Nv_Ne
%                                          min_Nv/Ne - min_Nv_Ne
%                         max_#_of_nodes_on_one_edge - max_on_edge
%                              intialization_inf_eps - eps
%                                             scaler - scaler
%                                              step# - stepN
%                                    iteration_times - Nit

% clear all;
%% set const
% for generating graphs
max_Nv = 20;
min_Nv = 10;
max_Nv_Ne = 5;
min_Nv_Ne = 3;
max_on_edge = min(max_Nv, ceil ( max_Nv_Ne * 1.5 ));

% for initialization
eps = 0.001;

% for simulation
Nit = 50;

% others
K = 100;
% global alp;
% alp = 2;
% scaler = 10;
% stepN = 1000;

%% generate hypergraph
randomG = 1;
unchanged = 1;

if unchanged == 0
    if randomG
        
        Nv = randi([min_Nv,max_Nv]);
        Ne = randi([ceil(Nv/max_Nv_Ne),ceil(Nv/min_Nv_Ne)]);
        G = zeros(Ne,Nv);
        for i = 1:Nv
            G(randi(Ne),i) = 1; % every node should belong to at least one edge
        end
        
        for itt = 1:Ne
            n_on_j = min(Ne,randi(max_on_edge)); % how many nodes choosen by this edge j
            perm = randperm(Nv);
            on_j = perm(1:n_on_j);
            G(itt,on_j) = 1;
        end
        
    else
        % mannually generate
        %     Nv = 4;
        %     Ne = 4;
        %     G = zeros(4,4);
        %     G(1,1) = 1;
        %     G(1,2) = 1;
        %     G(2,2) = 1;
        %     G(2,3) = 1;
        %     G(3,3) = 1;
        %     G(3,4) = 1;
        %     G(4,4) = 1;
        %     G(4,1) = 1;
        
        
        Nv = 4;
        Ne = 3;
        G = zeros(Ne,Nv);
        G(1,1) = 1;
        G(1,2) = 1;
        G(1,3) = 1;
        G(2,1) = 1;
        G(2,4) = 1;
        G(3,2) = 1;
        G(3,4) = 1;
    end
    
end

%% initialization
randomI = 1;

if randomI
    C = eps * rand(Nv,1);
    
else
    % mannually generate
    C = 0.1 * ones(Nv,1);
    C(5) = 500;
end



%% evolve
% result def
C_nom = zeros(Ne,1);
P_ser = zeros(Ne,Nv);
P_suc = zeros(Nv,1);
U = zeros(Nv,1);
% history record def
Crecord = zeros(Nit,Nv);
P_suc_record = zeros(Nit,Nv);

% running
for it = 1:Nit
    [C_nom, P_ser, P_suc, U] =  computeU(Ne,Nv,G,C);
    Crecord(it,:) = C';
    P_suc_record(it,:) = P_suc;
    C = P_suc * K;
end


%% PLOT
% plot(Crecord(1:25,:));
% figure;
% plot(sum(Crecord(1:25,:),2),'r*-');












% Urecord = zeros(Nit,Nv);
% Crecord = zeros(Nit,Nv);
% P_fail_record = zeros(Nit,Nv);
%
% K = 31;
%
% func = @(pos,C,Nv,G,x) sum ( 1 ./ ( G * [C(1:pos-1);x;C(pos+1:end)] * ones(1,Nv) ).* G);
% firstD = @(pos,C,Nv,G,x,K) K - sum( func(pos,C,Nv,G,x) .* [zeros(1,pos-1),1,zeros(1,Nv-pos)]);
%
% for it = 1:Nit
%
%     [C_nom, P_ser, P_suc, U] =  computeU(Ne,Nv,G,C);
%
%     Urecord(it,:) = U';
%     Crecord(it,:) = C';
%
%     P_fail = 1 - P_suc;
%     P_fail_record(it,:) = P_fail';
%
%     C_after = zeros(Nv,1);
%     for i = 1:Nv
%     %    XX(i) = sum( 1 / ( ( C_nom - C(i) * ones(Ne,1)) .* G(:,i) ));
%         %XX = sum( C_nom - C(i) * ones(Ne,1)) .* G(:,i));
%         %XX(XX==0) =
%         if K > sum( 1 / ( ( C_nom - C(i) * ones(Ne,1)) .* G(:,i) ))
%             C_after(i) = 0;
%         else
%             C_after(i) = binarySearch( @(x) firstD(i,C,Nv,G,x,K),[eps,10000]);
%         end
%     end
%     C = C_after;
% end
%
% % %% round by round
% % C_nom = zeros(Ne,1);
% % P_ser = zeros(Ne,Nv);
% % P_suc = zeros(Nv,1);
% % U = zeros(Nv,1);
% %
% % Urecord = zeros(Nit,Nv);
% % Crecord = zeros(Nit,Nv);
% % P_sucrecord = zeros(Nit,Nv);
% %
% % for it = 1:Nit
% %
% %     %     C_nom = G * C;
% %     %
% %     %     C_nom(C_nom==0) = -1;
% %     %     extended_C_nom = C_nom * ones(1,Nv);
% %     %     P_ser = G .* ( ones(Ne,1) * C'./ extended_C_nom );
% %     %
% %     %     p_temp = P_ser;
% %     %     p_temp(p_temp==0) = 1;
% %     %     P_suc = (1 - prod(p_temp))';
% %     %
% %     %     for i = 1:Nv
% %     %         diff = 0;
% %     %         for j = 1:Nv
% %     %             for k = 1:Ne
% %     %                 if G(k,i) == 1 && G(k,j) == 1
% %     %                     diff = diff + (P_suc(j) - P_suc(i)) ^ 2;
% %     %                     break;
% %     %                 end
% %     %             end
% %     %         end
% %     %         U(i) = 1/(alp + diff);
% %     %     end
% %
% %     [C_nom, P_ser, P_suc, U] =  computeU(Ne,Nv,G,C);
% %
% %     Urecord(it,:) = U';
% %     Crecord(it,:) = C';
% %     P_sucrecord(it,:) = P_suc';
% %
% %     C_after = zeros(Nv,1);
% %
% % %     for i = 1:Nv
% % %         cc = C;
% % %         maxU = U(i);
% % %         bestC = C(i);
% % %         for itt = 1:1000
% % %             cc(i) = 0.001*i;
% % %             [C_nom, P_ser, P_suc, uu] =  computeU(Ne,Nv,G,cc);
% % %             if uu > maxU
% % %                 bestC = cc(i);
% % %                 maxU = uu;
% % %             end
% % %         end
% % %         C_after(i) = bestC;
% % %     end
% %
% %     C_after = P_suc * alp;
% %
% %     C = C_after;
% % end
