%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com, https://github.com/lucasrm25)
%   - 
%   -

%   Generate fancy plots (fp)
%------------------------------------------------------------------

classdef fp
    properties
    end
    
    methods(Static)
        
        function savefig(fname,varargin)
            
            p = inputParser;
          
            defaultFormat  = 'epsc';
            expectedFormat = {'epsc','svg','png','jpg'};
            addParameter(p,'format',defaultFormat, @(x) any(validatestring(x,expectedFormat)));
            addParameter(p,'fighandle',gcf, @(x) isa(x,'matlab.ui.Figure'));
            addParameter(p,'folder',fullfile(pwd,'images'), @(x) isa(x,'char'));      
            parse(p,varargin{:});
            
            if ~ exist(p.Results.folder,'dir'), mkdir(p.Results.folder); end
            % set(p.Results.fighandle.CurrentAxes,'LooseInset',p.Results.fighandle.CurrentAxes.TightInset)
            saveas(p.Results.fighandle, fullfile(p.Results.folder,fname), p.Results.format)
        end
        
        function fig = f()
            fig = figure('Color','white');
            hold on, grid on;
        end
        
        function latexcode = m2latex(matrix)
            if ~isa(matrix,'sym')
                matrix = sym(matrix);
            end
            latexcode = latex(vpa(simplify(matrix)));
            if numel(latexcode)>=6 && strcmp(latexcode(6),'(') && strcmp(latexcode(end),')')
                latexcode(6) = '[';
                latexcode(end) = ']';
            end
            clipboard('copy',latexcode);
        end
        
        function str = figpos()
            a = gcf;
            pos = a.Position;
            str = ['figure(''Color'',''white'',''Position'',[' num2str(pos) ']);'];
            clipboard('copy',str);
        end
        
        % varargin{1}: line opacity
        function cl = getColor(n, varargin)
            colrs = lines(max(n));
            if nargin >= 2
                opacity = varargin{1};
            else
                opacity = [];
            end
            cl = [colrs(n,:), repmat(opacity,numel(n),1)];
        end
    end
end
