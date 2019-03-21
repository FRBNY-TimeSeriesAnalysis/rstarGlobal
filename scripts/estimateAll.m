%% Estimates specifications in paper. To estimate the appendix figures, switch 
% estimateAppendices = 1. Results are saved in the results folder.

estimateAppendices = 0;

% Specifications used in main body of paper
MainModel1              % Baseline model
MainModel1_var01        % Baseline model with disperse prior of variance of innovation to trend
MainModel2              % Convenience yield model
MainModel3              % Consumption model


if estimateAppendices == 1
    % Alterantive priors for variance of innovation to trend
    MainModel1_var02  
    MainModel1_var05
    MainModel1_var10
    MainModel1_var25
    MainModel1_var50
    
    MainModel1_unrestr  % Unrestricted loading for real rate
    MainModel1_ReR      % Baseline model with real exchange rate
    MainModel1_df50     % Baseline model with df = 50
    MainModel1_1950     % Baseline model, estimation starting 1950
end