function MetaViewer(varargin)
%MetaViewer Tree style metadata viewer
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Uses the undocumented matlab uitree control to interactively display
% image metadata.
%
% INPUTS:
%   filename/meta structure  -  optional : filename or meta structure
%                               missing  : user is prompted for filename
%
% OUTPUTS:
%   None
%
% VERSION:
%   1.0
%     - Tim Cox 20110523
%     - initial version, liberally borrowed some code off the internet
%   2.0
%     - Wade Schwartzkopf 20140121
%     - Changes introduced in the release of MATLAB 2014a broke the
%       previous version of this function.  Rewrote the same functionality
%       to work on more versions of MATLAB (tested on 2011a-2014a).

%determine if input is filename or structure
if ~isempty(varargin)&&isstruct(varargin{1})
    meta = varargin{1};
else
    if ~isempty(varargin) % Assume filename was passed
        filename = varargin{1};
    else % Otherwise ask
        %launch file browser
        %load last path
        if ispref('matlab_sar_toolbox','last_used_directory')
            pathstr = getpref('matlab_sar_toolbox','last_used_directory');
            if ~ischar(pathstr)||~exist(pathstr,'dir')
                pathstr = pwd;
            end
        else
            pathstr = pwd;
        end
        
        %get filename
        [filename, pathstr] = uigetfile( sar_file_extensions( {'complex','phd'} ), 'Open File',pathstr);
        if isnumeric(filename)
            return;
        end
        
        filename = strcat(pathstr,filename);
        
        setpref('matlab_sar_toolbox','last_used_directory',pathstr); %store path
    end
    %load meta
    try
        reader_obj = open_ph_reader(filename); % Phase history
    catch
        reader_obj = open_reader(filename); % Complex image
    end
    if iscell(reader_obj)
        reader_obj = reader_obj{1};
    end
    meta = reader_obj.get_meta();
    reader_obj.close();
end

% Draw figure
h = figure;
set(h,'Name','MetaViewer','MenuBar','none');

% Build tree
root = uitreenode('v0', 'string_id_we_dont_use', '', [], false); % Root node is empty dummy node
set(root,'UserData',meta); % Metadata is passed through uitreenodes
current_meta = struct(); % Variable visible by all callbacks
% If we wanted to build whole tree first:
% build_full_tree(root);
% tree = uitree('v0', h, 'Root', root);
% But we will use lazy loading:
tree = uitree('v0', h, 'Root', root, 'ExpandFcn', @getChildNodes); % Other nodes will only be added as required for expansion
set(tree, 'NodeWillExpandCallback', @getNodeUserData); % Called before each expansion
drawnow; % Must do this first for position setting to work on next line
set(tree, 'Units', 'normalized', 'position', [0 0 1 1]); % Whole figure

tree.expand(root); % Start with first layer of hierarchy exposed

%% Lazy loading
% We could create the entire tree structure up front traversing the entire
% metadata structure.  However, large, deep structures can take a long time
% to traverse completely, when the user likely only wants to see a handful
% of fields.  Therefore, we use "lazy loading"; that is, we only add nodes
% to the tree when expanding requires them.
% It would make most sense to put all of the "lazy loading" code in the
% ExpandFcn function, but we don't know how to get the handle to the
% uitreenode object being expanded within the ExpandFcn function.  The
% NodeWillExpandCallback function (which is called each time before a node
% expansion) does have access to the node handle, so we use it to pass data
% to the ExpandFcn (called only the first time a node is expanded to
% determine its children) through a shared variable.
    function node_to_expand = getNodeUserData(tree,event) % NodeWillExpandCallback
        node_to_expand = event.getCurrentNode;
        current_meta = get(node_to_expand,'UserData');
        if ~isempty(current_meta)
            % Metadata will be split up and passed to children in ExpandFcn
            % which will be called immediately after this, so no need to keep.
            set(node_to_expand,'UserData',[]);
        end
    end

    function child_nodes = getChildNodes(tree,value) % ExpandFcn callback
        fnames = fieldnames(current_meta); % Assumes current_meta was already set by getNodeUserData
        pth = [fileparts(mfilename('fullpath')) filesep];
        count = 1;
        
        for i = 1:numel(fnames)
            val = current_meta.(fnames{i});
            if ~isstruct(val) && ~iscell(val) % Leaf node
                [description, iconpath] = handle_leaf(fnames{i}, val);
                iconpath = [pth,iconpath];
                child_nodes(count) = uitreenode('v0',fnames{i},description,iconpath,true);
                count = count + 1;
            elseif isscalar(val) % Not a leaf node.  Recursively process farther.
                iconpath = [pth,'struct_icon.GIF'];
                child_nodes(count) = uitreenode('v0',fnames{i},fnames{i},iconpath,false);
                count = count + 1;
                set(child_nodes(end),'UserData',val);
            else % Array of nodes
                for j = 1:numel(val)
                    if iscell(val) % Is there a cleaner way to do this?
                        val_i = val{j};
                    else
                        val_i = val(j);
                    end
                    fnames_disp = strcat(fnames{i}, '(', num2str(j),')');
                    if isstruct(val_i) % Node has children
                        iconpath = [pth,'structarray_icon.GIF'];
                        child_nodes(count) =  uitreenode('v0', fnames_disp, fnames_disp, iconpath, 0);
                        count = count + 1;
                        set(child_nodes(end),'UserData',val_i);
                    else % Leaf node (we assume no arrays of arrays will be found)
                        [description, iconpath] = handle_leaf(fnames{i}, val_i);
                        iconpath = [pth,iconpath];
                        child_nodes(count) = uitreenode('v0',fnames_disp,description,iconpath,true);
                        count = count + 1;                        
                    end
                end
            end
        end
        
        if count==1
            child_nodes=[];
        end
    end

    % Helper function handles leaf node processing
    function [description, iconpath] = handle_leaf(fname, val)
        if isnumeric(val)
            iconpath = 'double_icon.GIF';
            if any(strcmpi(fname,{'CollectStart','DateTime'}))
                % Special handling for date numbers
                array_str = datestr(val);
            else
                % Handle (2-D) arrays of numbers
                array_str = num2str(val).';
                array_str(end+1,:) = ';';
                array_str(end) = '';
            end
            description = sprintf('%s: %s',fname,array_str);
        elseif ischar(val)
            iconpath = 'char_icon.GIF';
            description = sprintf('%s: %s',fname,val);
        elseif islogical(val)
            iconpath = 'logic_icon.GIF';
            if val
                description = sprintf('%s: true',fname);
            else
                description = sprintf('%s: false',fname);
            end
        elseif isobject(val)
            iconpath = 'obj_icon.GIF';
            description = sprintf('%s: unknown',fname);
        else
            iconpath = 'unknown_icon.GIF';
            description = sprintf('%s: unknown',fname);
        end
    end

    % If we wanted to build whole tree at once, we could do it like this:
    function build_full_tree(node)
        current_meta = get(node,'UserData');
        child_nodes = getChildNodes(); % Uses current_meta
        for i = 1:numel(child_nodes)
            if ~isempty(get(child_nodes(i),'UserData'))
                build_full_tree(child_nodes(i));
            end
            node.add(child_nodes(i));
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////