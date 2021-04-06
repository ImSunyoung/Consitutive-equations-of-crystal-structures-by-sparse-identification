%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%           Discovering constiutive models for hexagonal crystals
% 
%                                                                2021-04-06
%                                                      imsunyoung@snu.ac.kr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                      
clc; clear all; close all;
%% Load data
addpath('data')
type = 1;
materials = {'GaN', 'ZnO'};
refer_data = sprintf('Reference_data_hexa_%s.mat',materials{type});
load(refer_data);

%%%%%%%%%%%%%%%%%%%%%%% reference data information %%%%%%%%%%%%%%%%%%%%%%%% 
% column 1 : name of loading path (total # of loading path : 20)
% column 2 : Green-Lagrangian strain [E11 E22 E33 E23 E13 E31]
% column 3 : second Piola-Kirchhoff stress [S11 S22 S33 S23 S13 S31]
% column 4 : tangent modulus matrix [C11 C12 C13 C14 C15 C16 C22 C23...C66]
% column 5 : deformation gradient [F11 F12 F13 F21 F22 F23 F31 F32 F33]

%% Discovering constitutive model by sparse identification 
%%% Construct feature_library(theta);
% tangent_modulus_componenets_set(y_refer);
% [Ci; Cij; Cijk] => [i; ij; ijk](material_ceoff_candidate)  
train_loading_set = [1:19]; % among 20 loading path
[ theta, y_refer, material_coeff_candidate ] = func_feature_library(data, train_loading_set);

%%% <<< Least square >>> %%%
coeff_least_sq = theta \ y_refer;

full_coeff_least_sq = material_coeff_candidate;
full_coeff_least_sq(:,4) = coeff_least_sq;
short_coeff_least_sq = full_coeff_least_sq(full_coeff_least_sq(:,4)~=0,:);

% Least square error compute 
y_pred = theta*coeff_least_sq;
error = abs(y_refer-y_pred)./abs(y_refer);
avg_error_least_sq = sum(error)/length(y_refer)*100;


%%% <<< Sparse identification with Energy criteria >>> %%%
lambda1 = 0.05;
lambda2 = 0.005;
[coeff_sparse_id] = func_sparse_id_energy_cutoff(theta, y_refer, lambda1, lambda2);

full_coeff_sparse_id = material_coeff_candidate;
full_coeff_sparse_id(:,4) = coeff_sparse_id;
short_coeff_sparse_id = full_coeff_sparse_id(full_coeff_sparse_id(:,4)~=0,:);

% Sparse identification error compute
y_pred = theta*coeff_sparse_id;
error = abs(y_refer-y_pred)./abs(y_refer);
avg_error_sparse_id = sum(error)/length(y_refer)*100;

% save_name = sprintf('data/constitutive_coeff_hexa_%s.mat', materials{type});
% save(save_name, 'short_coeff_sparse_id')
%% Predict second Piola-Kirchhoff stress and tangent modulus
% fid = 1: predict second PK stress
% fid = 2: predict second PK stress and tangent modulus  
fid = 2; 

predict_path = 1; % among 20 loading path
nStep = size(data{predict_path,2},1); 

%%% fitting results
fitting_PK2 = zeros(3,3,nStep);
fitting_tmd = zeros(6,6,nStep);

iter = 0;
for i = 1:nStep
    iter = iter + 1;
    matF = reshape(data{predict_path,5}(i,:),3,[]);
    [var_out] = func_hexagonal_predict(fid, short_coeff_sparse_id, matF);
    if fid == 1
        fitting_PK2(:,:,iter) = var_out{1};
    else
        fitting_PK2(:,:,iter) = var_out{1};
        fitting_tmd(:,:,iter) = var_out{2};
    end
end

%%% plot predicted second PK stress
plot_pred_stress

%%% plot predicted tangent modulus
if fid == 2
    plot_pred_tangent_modulus
end

    
    