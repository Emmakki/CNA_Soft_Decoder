%% SOFT_DECODER_GROUPE3.m
% =========================================================================
% Author: Bibard Grau Emma; Beny Guillaume; Bek Nada; Elmoujaouid Boutaina; 
% Date: 2023, Octobre
% =========================================================================

function c_cor = SOFT_DECODER_GROUPE3(c, H, p, MAX_ITER)
% SOFT_DECODER_GROUPE3 Fonction de décodage LDPC: Soft decision
% 
% c : un vecteur colonne binaire de dimension [N, 1]
% H : une matrice de parité de dimension [M, N] constituée de true et false
% p : un vecteur colonne de dimension [N, 1]
% MAX_ITER : un entier strictement positif
% 
% c_cor le vecteur colonne binaire de dimension [N, 1]

%% FAIRE LES VERIFICATION DES ENTREES

%1ere Etape
% -------------------------------------------------------------------------
% c -> f: Qij ou Qij(1) = Pi et q(0) = 1-Pi
[M, N] = size(H)
% Construire une matrice q_0, q_1, r_0 et r_1 de dimenssion [M, N]
fprintf('Etape 1')
q_1 = p.*H
q_0 = (1-q_1).*H

for inter = 1:MAX_ITER
    % 2e Etape
    % ---------------------------------------------------------------------
    % Calcul des Rij(0) et Rij(1) de chaque f en fonction des q reçus
    fprintf('Etape 2')
    for i = 1:M
        for j = 1:N
            r_0(i,j) = check_node_message(j,q_1,H);
            r_1(i,j) = 1-r_0(i,j);
        end
    end
    r_0
    r_1
    
    % 3e Etape
    % ---------------------------------------------------------------------
    % Mise à jour des Qij en fonction des Qij reçus
    % /!\ Kij est choisi de tel sorte que Qij(0)+Qij(1) = 1
    fprintf('Etape 2')
    for i = 1:M
        for j = 1:N
            q_0(i,j), q_1(i,j) = variable_node_message (i,r_0, r_1, H, p(i));
        end
    end
    q_0
    q_1
    
    % Mise à jour de l'estimation avec Qi
    
    c_cor = estimation(r_0,r_1, p);
    
    % Test parity check sinon étape 2 sur MAX_ITER
    % wc << n et wr << m
    %% A COMPRENDRE
    parity_check = equation_parity_check();
    if parity_check == 1
        return
    end
    
end

end

function c_node_res = check_node_message (j,q_1,H)
% Calcule le message réponse r du neud f_j

prod = 1;
for i = 1:size(H,2)
    if H(i,j)==1 && i ~= j
        prod = prod*(1-2*q_1(i,j));
    end
c_node_res = (1/2) + (1/2)*prod;
end

end

%% A RETRAVAILLER
function [v_node_res0, v_node_res1] = variable_node_message(i, r_0, r_1, H, p)
% Calcule le message réponse q du neud c_i

prod0 = 1;
prod1 = 1;

for j = 1:size(H,1)
    if H(i,j)==1 && i ~= j
        prod0 = prod0*r_0(i,j);
        prod1 = prod1*r_0(i,j);
    end
end

v_node_res0 = (1-p)*prod0;
v_node_res1 = p*prod1;

K = 1/(v_node_res1+v_node_res0);

v_node_res0 = K*v_node_res0;
v_node_res1 = K*v_node_res1;

end

%% A RETRAVAILLER
function c_est = estimation (r_0,r_1, p)
% Détermine une estimation du mot code c
c_est = zeros(length(p),1);
for i = 1:length(p)
x0 = 1;
x1 = 1;
    for j = 1:size(r_0,1)
        x0 = x0*r_0(i,j);
        x1 = x1*r_1(i,j);
    end
qi_0 = (1-p(i))*x0;
qi_1 = p(i)*x1;

K = 1/(qi_0+qi_1);

qi_0 = K*qi_0;
qi_1 = K*qi_1;

if qi_1>qi_0
    c_est(i) = 1;
else
    c_est(i) = 0;

end
end
end

%% A CODER
function parity_check = equation_parity_check ()
parity_check = 1;
end