function [] = changelayout(varargin)
% LYT Quickly load or save a MATLAB desktop layout from the command line.
%   LYT by itself changes the MATLAB desktop layout to the default layout,
%   nominally 'Default'. This default layout can be changed in the code.
% 
%   LYT NAME changes the layout to the saved layout named NAME.
% 
%   LYT NAME -SAVE saves the current layout as NAME. If NAME is already a
%   saved layout, it will be overwritten without warning!
% 
%   For multi-word NAME, use LYT(NAME) syntax.
% 
%   See also PREFDIR, PREFERENCES, 
%     http://undocumentedmatlab.com/blog/setting-the-matlab-desktop-layout-programmatically/.

%   Copyright 2013 Sky Sartorius
%   Author contact: mathworks.com/matlabcentral/fileexchange/authors/101715



%% %%%%%%%%%%%%%%% User-definable default desktop layout %%%%%%%%%%%%%%%%%%
defaultLayoutName = 'Default';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Process inputs
saveStringIndex = strcmpi('-save',varargin);
if any(saveStringIndex)
    action = 'saveLayout';
    varargin = varargin(~saveStringIndex);
    % Don't care if there were multiple -save strings provided.
    if length(varargin) ~= 1
        error('LYT requires exactly one layout name for saving.')
    end
    layoutName = varargin{1};
    % NB: no use of default name allowed for saving.
else
    switch length(varargin);
        case 0
            layoutName = defaultLayoutName;
        case 1
            layoutName = varargin{1};
        otherwise % > 1
            error('LYT must be provided with maximum one layout name.')
    end
    action = 'restoreLayout';
end

%% Get desktop. Do action.
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.(action)(layoutName);


% Revision history
%{
2013-10-22 Added -save option and restructure whole function w/ input
checks
%}