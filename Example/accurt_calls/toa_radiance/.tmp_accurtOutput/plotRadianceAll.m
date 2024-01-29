% This Matlab script plots the AccuRT computed angular radiance
% distribution at a given run for ALL wavelengths and depths.
% WARNING: while this provides a good overview of the output,
% it also may produce many plots depending on how many wavelengths
% and depths are specified (wavelengths*depths number of plots).

clear all;

% Begin input

fileName     = 'radiance.txt';
runNo        = 1;

% End input

data = readRadiance(fileName);

for depthNo=1:data(runNo).nDepths
    for wavelengthNo=1:data(runNo).nWavelengths
    
    figure;

    x = data(runNo).polarAngles;
    y = data(runNo).azimuthAngles;
    z = squeeze(data(runNo).radiance(depthNo,wavelengthNo,:,:))';

    xlab = 'Polar angle [degree]';
    ylab = 'Azimuth angle [degree]';
    if length(y)==1
        tmp = x;
        x = y;
        y = tmp;
        tmp = xlab;
        xlab = ylab;
        ylab = tmp;
    end

    h = bar3(z);
    shading interp
    for k = 1:length(h)
        zdata = get(h(k),'ZData');
        set(h(k),'CData',zdata)
        set(h,'EdgeColor','k') 
    end
    hcb = colorbar;
    set(get(hcb,'title'),'string','[W m^{-2} nm^{-1} sr^{-1}]')
    rx = 1:ceil(length(x)/8):length(x);
    ry = 1:ceil(length(y)/8):length(y);
    xl = round(x(rx)*10)/10;
    yl = round(y(ry)*10)/10;  
    set(gca,'xtick',rx)
    set(gca,'xticklabel',xl)
    set(gca,'ytick',ry)
    set(gca,'yticklabel',yl)
    title(['run = ',num2str(runNo),', depth = ', num2str(data(runNo).depths(depthNo)),' m', ...
       ', wavelength = ',num2str(data(runNo).wavelengths(wavelengthNo)),' nm', ...
       ', streams = ', num2str(data(runNo).nStreams)])
    xlabel(xlab,'fontname','times')
    ylabel(ylab,'fontname','times')
    zlabel('Radiance [W m^{-2} nm^{-1} sr^{-1}]','fontname','times')
    set(gca,'fontname','times')

    end
end
