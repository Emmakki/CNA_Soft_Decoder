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

% Réccupération de la taille de la matrice H
[M, N] = size(H);

% Initialisation des matrices contenant les messages transmits entre les
% noeuds (v_nodes et c_nodes)
v_nodes_msg_0 = -1*ones(M,N);
v_nodes_msg_1 = -1*ones(M,N);

% Envoie du premier message 
% Des variable nodes (v_nodes) vers les check nodes (c_nodes)
% Avec les probabilités p pour c(i) == 1 sachant le signal reçu
for j = 1:M
    for i = 1:N
        if H(j,i) == 1
            v_nodes_msg_0(j,i) = 1-p(j);
            v_nodes_msg_1(j,i) = p(j);
        end
    end
end

% 2e Etape
% ----------------------------------------------------------------------
for inter = 1:MAX_ITER
    % Calcul des messages réponse des c_nodes vers les v_nodes
    c_nodes_msg_0 = reponse_check_nodes (H,v_nodes_msg_1);
    c_nodes_msg_1 = 1-c_nodes_msg_0;
    
    % 3e Etape
    % ------------------------------------------------------------------
    
    % Calcul des messages réponse des v_nodes vers les c_nodes
    [v_nodes_msg_0, v_nodes_msg_1] = reponse_variable_nodes(H,p,c_nodes_msg_0, c_nodes_msg_1);

    % Mise à jour de l'estimation avec Qi
    c_cor = estimation(c_nodes_msg_0,c_nodes_msg_1, p, H);
    
    % Test du parity check sinon retour à l'étape 2 sur MAX_ITER
    parity_check = equation_parity_check(H, c_cor);
    if parity_check == 1
        return
    end
    
end
end

function c_nodes_rsp = reponse_check_nodes (H,v_nodes_msg_1)
% Calcul la matrice contenant les messages des c-nodes (lignes) vers les
% v-nodes (colonnes)
[M, N] = size(H);

c_nodes_rsp = -1*ones(M,N);

for j = 1:M
    for i = 1:N
        if H(j,i) == 1
            prod = 1;
            for x = 1:size(H,2)
                if H(j,x) == 1 && x ~= i
                    prod = prod*(1-2*v_nodes_msg_1(j,x));
                end
            end
            c_nodes_rsp(j,i) = (1/2) + (1/2)*prod;
        end
    end
end
end

function [v_nodes_rsp_0, v_nodes_rsp_1] = reponse_variable_nodes(H,p,c_nodes_msg_0, c_nodes_msg_1)
%Initialisation des matrices
[M, N] = size(H);
v_nodes_rsp_0 = -1*ones(M,N);
v_nodes_rsp_1 = -1*ones(M,N);

%Calcule des réponses des v-nodes i aux c-nodes j
for j = 1:M
    for i = 1:N
        if H(j,i) == 1

            % Produit des messages des c-nodes reçus par le v-nodes i
            % A l'exception du messages du c-node j auquel on envoie
            c_prod_0 = 1;
            c_prod_1 = 1;
            for x = 1:size(H,1)
                if H(x,i) == 1 && x ~= j
                    c_prod_0 = c_prod_0*c_nodes_msg_0(x,i);
                    c_prod_1 = c_prod_1*c_nodes_msg_1(x,i);
                end
            end
            v_nodes_rsp_0(j,i) = (1-p(j))*c_prod_0;
            v_nodes_rsp_1(j,i) = p(j)*c_prod_1;

            %Calcul de la constante K
            % Où v_nodes_rsp_0(j,i) + v_nodes_rsp_1(j,i) = 1
            K = 1/(v_nodes_rsp_0(j,i)+v_nodes_rsp_1(j,i));

            % Valeur finale de la réponse du v-node i au c-node j
            v_nodes_rsp_0(j,i) = K*v_nodes_rsp_0(j,i);
            v_nodes_rsp_1(j,i) = K*v_nodes_rsp_1(j,i);

        end
    end
end
end

function c_est = estimation (c_nodes_msg_0,c_nodes_msg_1, p, H)
% Détermine une estimation du mot code c

%Initialisation de la matrice colonne corrigeant le mot code
[M, N] = size(H);
c_est = zeros(length(p),1);

% Calcule de Q avec l'aide chaque colonne des matrices messages des c-nodes
for i = 1:N
q_prod_0 = 1;
q_prod_1 = 1;
    for j = 1:M
        if H(j,i) == 1
            q_prod_0 = q_prod_0*c_nodes_msg_0(j,i);
            q_prod_1 = q_prod_1*c_nodes_msg_1(j,i);
        end
    end
q_0 = (1-p(i))*q_prod_0;
q_1 = p(i)*q_prod_1;

%Calcul de la constante K
K = 1/(q_0+q_1);

% Valeur finale pour Q associé à la composante j du mot code
q_0 = K*q_0;
q_1 = K*q_1;

% Correction du mot code 
if q_1>q_0
    c_est(i) = 1;
else
    c_est(i) = 0;
end
end
end

%% A CODER
function parity_check = equation_parity_check (H, c)
% Vérifie si le mot code c vérifie les équations de parité de la matrice H
vecteur = mod(H*c,2);
if sum(vecteur) == 0
    parity_check = 1;
else
    parity_check = 0;
end
end