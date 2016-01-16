function [C_nom, P_ser, P_suc, U] =  computeU(Ne,Nv,G,C)
% global alp;

% C_nom = zeros(Ne,1);
% P_ser = zeros(Ne,Nv);
% P_suc = zeros(Nv,1);
U = zeros(Nv,1);

%% compute C_nom and possibilities
C_nom = G * C;
C_nom_temp = C_nom;
C_nom_temp(C_nom==0) = -1;
extended_C_nom = C_nom_temp * ones(1,Nv);
P_ser = G .* ( ones(Ne,1) * C'./ extended_C_nom );

p_temp = P_ser;
p_temp(p_temp==0) = 1;
P_suc = (1 - prod(p_temp))';

%% compute U
% total diff between this node and neighbors
% for i = 1:Nv
%     diff = 0;
%     for j = 1:Nv
%         for k = 1:Ne
%             if G(k,i) == 1 && G(k,j) == 1
%                 diff = diff + (P_suc(j) - P_suc(i)) ^ 2;
%                 break;
%             end
%         end
%     end
%     U(i) = 1/(alp + diff);
% end

% new utility
% U = - (log(1-P_suc) + 31 * C);


end