function [ground_track_fig] = GroundTrackPlotSetup()

    ground_track_fig = figure();
    hold on
    imData = imread('no_clouds_4k.jpg'); 
    xlim([-180,180])
    ylim([-90,90])
    h = image(xlim,-ylim,imData); 
    uistack(h,'bottom')
    xlabel('Longitude °')
    ylabel('Latitude°')
    title('Ground Trace')
    ax = gca;
    ax.FontSize = 20;

end