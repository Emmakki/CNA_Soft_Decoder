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

% FAIRE LES VERIFICATION DES ENTREES

%1ere Etape
% -------------------------------------------------------------------------

% Initialisation du mot code c_cor qui va être corigé à chaque itération
c_cor = c;
%Réccupération de la taille des matrices
[M, N] = size(H);

% Envoie du premier message 
% Des variable nodes (v_nodes) vers les check nodes (c_nodes)
% Initialisation des matrices messages
v_nodes_msg_0 = ones(M,N);
v_nodes_msg_1 = ones(M,N);
c_nodes_msg_0 = ones(M,N);
c_nodes_msg_1 = ones(M,N);

% Remplissage de la matrice message des variable nodes
% Avec les probabilités p pour c(i) == 1 sachant le signal reçu
v_nodes_msg_0 = (1-p).*H;
v_nodes_msg_1 = p.*H;

% 2e Etape
% -------------------------------------------------------------------------
for inter = 1:MAX_ITER
    % Calcul des messages réponse des c_nodes vers les v_nodes
    c_nodes_msg_0 = reponse_check_nodes (H,v_nodes_msg_1);
    c_nodes_msg_1 = 1-c_nodes_msg_0;
    
    % 3e Etape
    % ---------------------------------------------------------------------
    % Calcul des messages réponse des v_nodes vers les c_nodes
    % /!\ Kij est choisi de tel sorte que Qij(0)+Qij(1) = 1
    for i = 1:M
        for j = 1:N
            q_0(i,j), q_1(i,j) = variable_node_message (i,r_0, r_1, H, p(i));
        end
    end

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

function c_nodes_rsp = reponse_check_nodes (H,v_nodes_msg_1)
c_nodes_rsp = zeros(size(H,1),size(H:2));
for j = 1:size(H,1)
    for i = 1:size(H:2)
        if H(j,i) == 1
            prod = 1;
            for p = 1:size(H,2)
                if H(j,p) == 1 && p ~= i
                    prod = prod*(1-2*v_nodes_msg_1(j,p));
                end
            end
            c_nodes_rsp(j,i) = (1/2) + (1/2)*prod;
        end
    end
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