% frm = getFrame(fig,roi)
%
% Figure window capture function that does not capture other things on the
% screen. Functions identical to getframe(fig,rect).
%
% Use the zbuffer or OpenGL renderer for this to work best.
% Example:
%    fig=figure('Renderer','zbuffer');
%    fig=figure('Renderer','OpenGL');
%    fig=figure();
%    set(fig,'Renderer','zbuffer');
%
% Usage:
%    frm=getFrame(fig);
%    writeVideo(videoWriter,frm);
%
function frm = getFrame(fig,roi)
    if (nargin<1 || isempty(fig))
        fig=gcf();
    end
    if (nargin<2)
        roi=[];
    end
    if strcmpi(get(fig,'Type'),'axes')
        fig=get(fig,'Parent');
    end
    
    switch lower(get(fig,'Renderer'))
        case 'zbuffer'
            device='-Dzbuffer';
        case 'opengl'
            device='-DOpenGL';
        otherwise
            errorMessage=sprintf('Renderer %s not supported for frame grabbing with getFrame(), use figure(''Renderer'',''zbuffer''). Switching to zbuffer now.',get(fig,'Renderer'));
            logMessage(errorMessage);
%             origRenderer=get(fig,'Renderer');
            set(fig,'Renderer','zbuffer');
            device='-Dzbuffer';
    end
    % Need to have PaperPositionMode be auto 
    origMode = get(fig, 'PaperPositionMode');
    set(fig, 'PaperPositionMode', 'auto');

    try
      cdata = print(fig,'-RGBImage');
    catch Exc
      cdata = hardcopy(fig, device, '-r0'); % Undocumented, but only option in old Matlab versions
    end
    cdata = cdata(1:(2*floor(end/2)),1:(2*floor(end/2)),:);
    
    if (~isempty(roi))
        cdata=cdata(max(1,roi(1)):min(end,roi(1)+roi(3)),max(1,roi(2)):min(end,roi(2)+roi(4)),:);
    end

    % Restore figure to original state
    set(fig, 'PaperPositionMode', origMode);
    
    frm = im2frame(cdata);
end
