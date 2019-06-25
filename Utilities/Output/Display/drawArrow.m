%
% ar = drawArrow(startPos, endPos, thickness)
%
% Draws an arrow on the current axes
%
function ar = drawArrow(startPos, endPos, thickness, varargin)
  if nargin < 1 || isempty(startPos)
    startPos = [0 0];
  end
  if nargin < 2
    endPos = [1 1];
  end
  tangential = endPos - startPos;
  orthogonal = [-tangential(2) tangential(1)];
  
  arrowLength = norm(tangential);
  if nargin < 3 || isempty(thickness)
    thickness = 0.02 * arrowLength;
  end
  headWidth = 3 * thickness;
  headLength = 2 * headWidth;
  
  P = bsxfun(@plus, startPos(:),...
    bsxfun(@times, tangential(:) / norm(tangential), [0 arrowLength-headLength arrowLength-headLength arrowLength arrowLength-headLength arrowLength-headLength 0]) );
  P = P + bsxfun(@times, orthogonal(:) / norm(orthogonal), [thickness thickness headWidth 0 -headWidth -thickness -thickness] ./ 2);
  
  ar = patch('XData', P(1,:), 'YData', P(2,:), varargin{:});
end