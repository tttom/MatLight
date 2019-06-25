% [back_aperture_diameter, objective_front_focal_length] = calcBackApertureDiameter(magnification, numericalAperture, tubeLensFocalLength)
%
% Calculates the diameter of the back aperture of an infinity-corrected microscope objective
% from the values typically written (or deduced from) the barrel
%
% Input:
%    magnification: the number in front of the 'X'
%    numerical_aperture: the number in front of the '/...X'
%    tube_lens_focal_length [m]: the most common tube length lens is 200 mm,
%          used for Nikon and Leica. Olympus uses 180 mm, and Zeiss 165 mm.
%
% Output:
%    back_aperture_diameter [m]: the back aperture diameter in meters
%    objective_front_focal_length [m]: the front focal length of the objective lens
function [back_aperture_diameter, objective_front_focal_length] = calcBackApertureDiameter(magnification, numerical_aperture, tube_lens_focal_length)
    if nargin<3
        tube_lens_focal_length = 200e-3;
    end
    
    % Assume infinity-corrected objective
    objective_front_focal_length = tube_lens_focal_length ./ abs(magnification);
    
    % Assume that the sine condition holds
%     % Don't assume the small angle approximation
%     sine_theta_image_space = numerical_aperture / magnification;
%     tangent_theta_image_space = sine_theta_image_space ./ sqrt(1 - sine_theta_image_space.^2);
%     back_aperture_diameter = 2 * tangent_theta_image_space .* tube_lens_focal_length;
    % Assume the small angle approximation (difference typically less than 1%)
    back_aperture_diameter = 2 * numerical_aperture .* objective_front_focal_length;
end