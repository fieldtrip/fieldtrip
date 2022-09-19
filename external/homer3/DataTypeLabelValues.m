function  dataTypeLabelValues = DataTypeLabelValues()
dataTypeLabelValues = struct(...
    'dOD','dOD', ...   % Change in optical density
    'mua','mua', ...   % Absorption coefficient
    'musp','musp', ...   % Scattering coefficient
    'HbO','HbO', ...   % Oxygenated hemoglobin (oxyhemoglobin) concentration
    'HbR','HbR', ...   % Deoxygenated hemoglobin (deoxyhemoglobin) concentration
    'HbT','HbT', ...   % Total hemoglobin concentration
    'HRF_HbO','HRF HbO', ...   % Hemodynamic response function for oxyhemoglobin concentration
    'HRF_HbR','HRF HbR', ...   % Hemodynamic response function for deoxyhemoglobin concentration
    'HRF_HbT','HRF HbT', ...   % Hemodynamic response function for total hemoglobin concentration
    'H2O','H2O', ...   % Water content
    'Lipid','Lipid', ...   % Lipid concentration
    'BFi','BFi', ...   % Hemodynamic response function for blood flow index (BFi)
    'HRF_BFi','HRF BFi' ...   % Hemodynamic response function for blood flow index (BFi)
);