function A = rot2ang(T)
% Get angles from rotation matrix.
% Direct copy from MNE-Python function mne.transform.rotation_angles by  
% Alxandre Gramfort and Christian Brodbeck 
% https://github.com/mne-tools/mne-python/blob/maint/0.15/mne/transforms.py#L75-L129

    x = atan2(T(3,2),T(3,3));
    c2 = sqrt(T(1,1)^2 + T(2,1)^2);
    y = atan2(-T(3,1), c2);
    s1 = sin(x);
    c1 = cos(x);
    z = atan2(s1 * T(1,3) - c1 * T(1,2), c1 * T(2,2) - s1 * T(2,3));
    
    A = [x,y,z];
    
end

% Original Python code
%     x = np.arctan2(m[2, 1], m[2, 2])
%     c2 = np.sqrt(m[0, 0] ** 2 + m[1, 0] ** 2)
%     y = np.arctan2(-m[2, 0], c2)
%     s1 = np.sin(x)
%     c1 = np.cos(x)
%     z = np.arctan2(s1 * m[0, 2] - c1 * m[0, 1], c1 * m[1, 1] - s1 * m[1, 2])
%     return x, y, z

