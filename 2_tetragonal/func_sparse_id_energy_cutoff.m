function [ Xi ] = func_sparse_id_energy_cutoff(Theta, y_refer, lambda1, lambda2)

numData = size(Theta,1);

% initial guess: Least-squares
Xi = Theta \ y_refer;  
% value per each coefficient 
energy_0 = Theta.*Xi';
% relative energy(abs(value)) per each coefficient 
energy_1 = abs(energy_0'./y_refer')';
% lambda1 : relative energy cutoff value
energy_12 =  (energy_1 > lambda1);
energy_2 = sum(energy_12,1)';

% sequantial process
for k = 1:10
    % lambda2 : cutoff ratio
    smallinds = (energy_2 < numData * lambda2);
    Xi(smallinds) = 0;
    energy_0(:,smallinds) = 0;
    
    % biginds : remained index of big ratio
    biginds = ~smallinds(:);
    
    Xi(biginds) = Theta(:,biginds) \ y_refer(:);
    energy_0(:,biginds) = Theta(:,biginds).*Xi(biginds)';
    energy_1 = abs(energy_0'./y_refer')';
    energy_12 =  (energy_1 > lambda1);
    energy_2 = sum(energy_12,1)';
    
end
