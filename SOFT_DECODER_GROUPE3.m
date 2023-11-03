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
v_nodes_msg_0 = -1*ones(N,M);
v_nodes_msg_1 = -1*ones(N,M);

% Envoie du premier message 
% Des variable nodes (v_nodes) vers les check nodes (c_nodes)
% Avec les probabilités p pour c(i) == 1 sachant le signal reçu
for i = 1:N
    for j = 1:M
        if H(j,i) == 1
            v_nodes_msg_0(i,j) = 1-p(i);
            v_nodes_msg_1(i,j) = p(i);
        end
    end
end

% 2e Etape
% ----------------------------------------------------------------------
for inter = 1:MAX_ITER
    % Calcul des messages réponse des c_nodes vers les v_nodes
    [c_nodes_msg_0, c_nodes_msg_1] = reponse_check_nodes (H,v_nodes_msg_1);
    
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

function [c_nodes_rsp_0, c_nodes_rsp_1] = reponse_check_nodes (H,v_nodes_msg_1)
% Calcul la matrice contenant les messages des c-nodes vers les
% v-nodes
[M, N] = size(H);

c_nodes_rsp_0 = -1*ones(M,N);
c_nodes_rsp_1 = -1*ones(M,N);

for i = 1:M
    for j = 1:N
        if H(i,j) == 1

            % Récupération de tout les messages reçus pour le c_nodes i
            v_nodes_i = v_nodes_msg_1(1:end,i);
            % Retire les v_nodes n'étant pas connectés au c_nodes
            p = v_nodes_i(v_nodes_i ~= -1);
            % Produit des messages reçu par le c_nodes sauf celui envoyé
            % par le v_nodes duquel on calcule le message
            p = prod(1-2*p)/(1-2*v_nodes_msg_1(j,i));

            c_nodes_rsp_0(i,j) = (1/2) + (1/2)*p;
            c_nodes_rsp_1(i,j) = 1-c_nodes_rsp_0(i,j);
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
for i = 1:N
    for j = 1:M
        if H(j,i) == 1

            % Récupération de tout les messages reçus pour le v_nodes i
            c_nodes_i_0 = c_nodes_msg_0(1:end,i);
            c_nodes_i_1 = c_nodes_msg_1(1:end,i);
            % Retire les c_nodes n'étant pas connectés au v_nodes
            produit_0 = c_nodes_i_0(c_nodes_i_0 ~= -1);
            produit_1 = c_nodes_i_1(c_nodes_i_1 ~= -1);
            % Produit des messages reçu par le v_nodes sauf celui envoyé
            % par le c_nodes duquel on calcule le message
            produit_0 = prod(produit_0)/(c_nodes_msg_0(j,i));
            produit_1 = prod(produit_1)/(c_nodes_msg_1(j,i));

            v_nodes_rsp_0(i,j) = (1-p(i))*produit_0;
            v_nodes_rsp_1(i,j) = p(i)*produit_1;

            %Calcul de la constante K
            % Où v_nodes_rsp_0(j,i) + v_nodes_rsp_1(j,i) = 1
            K = 1/(v_nodes_rsp_0(i,j)+v_nodes_rsp_1(i,j));

            % Valeur finale de la réponse du v-node i au c-node j
            v_nodes_rsp_0(i,j) = K*v_nodes_rsp_0(i,j);
            v_nodes_rsp_1(i,j) = K*v_nodes_rsp_1(i,j);

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

function parity_check = equation_parity_check (H, c)
% Vérifie si le mot code c vérifie les équations de parité de la matrice H
vecteur = mod(H*c,2);
if sum(vecteur) == 0
    parity_check = 1;
else
    parity_check = 0;
end
end