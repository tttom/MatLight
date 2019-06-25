% outputFilePath = saveColorMap(colorMap, outputFilePath)
%
% Function to create a color bar image for inclusion with false color
% figures.
%
% colorMap: a Nx3 matrix with the values of red, green, and blue, in the respective columns. default: the current color map
% outputFilePath: the file name (and path) to save the png, fig, and svg files. Default: colorMap in the current folder.
%
function outputFilePath = saveColorMap(colorMap, outputFilePath)
    if nargin<1 || isempty(colorMap),
        colorMap = colormap();
    end
    if nargin<2,
        outputFilePath = fullfile(pwd(),'colorMap');
    end
    
    scalingFactor = 4;
    
    horizontalPx = 16;
    
    fig = figure('Renderer','zbuffer','Position',[50 50 [64 256]*scalingFactor],'Color','none');
    ax = axes();
    img = repmat(permute(colorMap,[1 3 2]),[1 horizontalPx 1]);
    image([1:size(img,2)]./size(img,2),([1:size(img,1)]-1)./(size(img,1)-1),img);
    lab = ylabel('I [a.u.]');
    set(ax,'LineWidth',scalingFactor,'YDir','normal','XTick',[],'YTick',-10:.5:10,'TickDir','in');
    set([ax, lab],'FontSize',9*scalingFactor,'FontWeight','bold');
    axis(ax,'tight');
    set(ax,'Position',[0.67 0.05 0.25 0.9])
    
    % Save images
    saveas(fig,[outputFilePath,'.fig']);
    saveWithTransparency(fig,[outputFilePath,'.png']);
    plot2svg([outputFilePath,'.svg'],fig,'png');
    
    close(fig);
end