% This Matlab script plots the AccuRT computed angular radiance
% distribution at a given run, depth, and wavelength in 2-dimensions.
% Each plot represents the result for a particular azimuthal angle.

clear all;

% Begin input

fileName     = 'radiance.txt';
runNo        = 1;
depthNo      = 1; 
wavelengthNo = 1;

% End input

data = readRadiance(fileName);
nStreams = data(runNo).nStreams;

for azimuthNo=1:data(runNo).nAzimuthAngles
    
    figure

    x = data(runNo).polarAngles;
    y = data(runNo).azimuthAngles(azimuthNo);
    z = squeeze(data(runNo).radiance(depthNo,wavelengthNo,:,azimuthNo))';

    xlab = 'Polar angle [degree]';
    plot(x,z)
    title(['run = ',num2str(runNo),', depth = ', num2str(data(runNo).depths(depthNo)),' m', ...
       ', wavelength = ',num2str(data(runNo).wavelengths(wavelengthNo)),' nm, azimuth = ',num2str(y),' degrees', ...
       ', streams = ', num2str(nStreams)])
    xlabel(xlab,'fontname','times')
    ylabel('Radiance [W m^{-2} nm^{-1} sr^{-1}]','fontname','times')
    set(gca,'fontname','times')

end
