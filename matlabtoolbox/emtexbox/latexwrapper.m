function outfile = latexwrapper(file, command, varargin)
% function filestruct = latextablewrap(file, command, varargin)
% commands: start/stop/ add [type] [filename]

%   Coded by  Elmar Mertens, em@elmarmertens.com

% TODO: replace calls to !start with winopen

narginchk(0, 2 + length(varargin));

if nargin < 1
   evalin('caller', 'latexwrapper(wrap)')
   return
end

if ~isstruct(file)
   if ischar(file)
      file.name = file;
   else
      error('file argument needs to be structure or char!')
   end
else
   if ~isfield(file, 'name')
      error('file structure needs to contain field "name"')
   end
end

if ~isfield(file, 'dir')
   file.dir = '.';
else
   if isempty(dir(file.dir))
      mkdir(file.dir)
   end
end

if ~isfield(file, 'pagestyle')
   file.pagestyle = 'headings';
end
if ~isfield(file, 'dvidriver')
   if isunix
      file.dvidriver = 'dvips'; % 2015, works better this way at the Board
   else
%       file.dvidriver = 'dvipdfm';
      file.dvidriver = 'dvips';
   end
end

filename = fullfile(file.dir, strcat(file.name, '.tex'));

if ~isfield(file, 'title')
   file.title = [];
   file.author = [];
elseif ~isfield(file, 'author')
   file.author = [];
end

if nargin < 2 || isempty(command)
   command = 'open';
end

switch lower(command)
   case {'start', 'setup'}

      [file.id, message] = fopen(filename, 'wt');
      if file.id < 0
         error('Could not open file: "%s" message: "%s"', filename, message)
      end
      touchfile = sprintf('%s.touch', fullfile(file.dir, file.name));
      if exist(touchfile, 'file')
         warning('em:msg', 'wrap %s seems to be in use already!', fullfile(file.dir, file.name))
         %          fprintf('(you have 60 sec to abort this program)\n')
         %          pause(60)
      end
      fclose(fopen(touchfile, 'wt'));
      %       if ispc
      %           fclose(fopen(touchfile, 'wt'));
      %       else
      %           system(sprintf('touch %s', touchfile));
      %       end

      fprintf(file.id, '%s\n', '\documentclass[11pt]{article}');
      fprintf(file.id, '%s\n', '\usepackage{exscale, xspace, setspace, ifthen}');
      fprintf(file.id, '%s\n', '\usepackage{amsmath, amsbsy, amsfonts, dcolumn, booktabs} ');
      fprintf(file.id, '%s\n', '\usepackage{afterpage,lastpage}');
      fprintf(file.id, '%s\n', '\usepackage{lscape}');

      fprintf(file.id, '%s\n', '\usepackage{ifpdf}');
      fprintf(file.id, '%s\n', '\ifpdf');
      fprintf(file.id, '%s\n', '  \usepackage{graphicx,rotating,hyperref}');
      fprintf(file.id, '%s\n', '\else');
      fprintf(file.id, '%s\n', '  \usepackage[dvips]{graphicx,rotating,hyperref}');
      fprintf(file.id, '%s\n', '  \DeclareGraphicsExtensions{.eps}');
      fprintf(file.id, '%s\n', '\fi');
      
      fprintf(file.id, '%s\n', '\usepackage[latin1]{inputenc}');
      fprintf(file.id, '%s\n', '\usepackage{listings}');
      fprintf(file.id, '%s\n', '\usepackage{a4wide}');
      fprintf(file.id, '%s\n', '\newcommand*{\titlecaveat}[1]{\texttt{-----------------------------\\ #1 \\-----------------------------}}');
      fprintf(file.id, '%s\n', '\newcounter{XYZs}');
      fprintf(file.id, '%s\n', '\newcommand*{\XYZ}[1][XYZ]{\fbox{\texttt{\textbf{#1}}}\xspace\addtocounter{XYZs}{1}}');
      fprintf(file.id, '%s\n', '\newcommand*{\marXYZ}[1][XYZ]{\XYZ[#1]\marginpar{\fbox{\texttt{CHECK}}}}');
      fprintf(file.id, '%s\n', '\newcommand*{\ccol}[1]{\multicolumn{1}{c}{#1}}');
      fprintf(file.id, '%s\n', '\newcommand*{\rcol}[1]{\multicolumn{1}{r}{#1}}');
      fprintf(file.id, '%s\n', '\newcommand*{\lcol}[1]{\multicolumn{1}{l}{#1}}');
      fprintf(file.id, '%s\n', '\newcolumntype{.}[1]{D{.}{.}{#1}}');
      fprintf(file.id, '%s\n', '\newcommand*{\ZstarText}{$^{\ast\ast\ast}$, $^{\ast\ast}$ and $^{\ast}$ denote significance at the 1\%, 5\% respectively 10\% level.\xspace}');
      fprintf(file.id, '%s\n', '\newcommand*{\FloatHead}{\hrulefill}');
      fprintf(file.id, '%s\n', '\newcommand*{\FloatBottom}{\FloatHead}');
      % fprintf(file.id, '%s\n', '\newcommand{\imagesandtables}{}');
      fprintf(file.id, '%s\n', '\newcommand*{\inputTAB}[2]{');
      fprintf(file.id, '%s\n', '  \begin{table}[t]');
      fprintf(file.id, '%s\n', '  \caption{#2}');
      % fprintf(file.id, '%s\n', '  \begin{center}');
      fprintf(file.id, '%s\n', '  \input{#1}');
      % fprintf(file.id, '%s\n', '  \end{center}');
      fprintf(file.id, '%s\n', '  \begin{center} \texttt{File: #1} \end{center}');
      fprintf(file.id, '%s\n', '  \end{table}');
      fprintf(file.id, '%s\n', '  \clearpage');
      fprintf(file.id, '%s\n', '}');

      fprintf(file.id, '%s\n', '\newcommand*{\inputTABHERE}[2]{');
      fprintf(file.id, '%s\n', '  \begin{table}[h]');
      fprintf(file.id, '%s\n', '  \caption{#2}');
      fprintf(file.id, '%s\n', '  \begin{center}');
      fprintf(file.id, '%s\n', '  \input{#1}');
      fprintf(file.id, '%s\n', '  \end{center}');
      fprintf(file.id, '%s\n', '  \begin{center} \texttt{File: #1} \end{center}');
      fprintf(file.id, '%s\n', '  \end{table}');
      %      fprintf(file.id, '%s\n', '  \clearpage');
      fprintf(file.id, '%s\n', '}');
     
      fprintf(file.id, '%s\n', '\newcommand*{\inputSIDETAB}[2]{');
      fprintf(file.id, '%s\n', '  \begin{sidewaystable}[t]');
      fprintf(file.id, '%s\n', '  \caption{#2}');
      %       fprintf(file.id, '%s\n', '  \begin{center}');
      fprintf(file.id, '%s\n', '  \input{#1}');
      %       fprintf(file.id, '%s\n', '  \end{center}');
      fprintf(file.id, '%s\n', '  \begin{center} \texttt{File: #1} \end{center}');
      fprintf(file.id, '%s\n', '  \end{sidewaystable}');
      fprintf(file.id, '%s\n', '  \clearpage');
      fprintf(file.id, '%s\n', '}');
      fprintf(file.id, '%s\n', '\newcommand*{\inputFIG}[3]{');
      fprintf(file.id, '%s\n', '  \begin{figure}[t]');
      fprintf(file.id, '%s\n', '  \caption{#2}');
      fprintf(file.id, '%s\n', '  \includegraphics[width=\textwidth]{#1}');
      fprintf(file.id, '%s\n', '  \begin{footnotesize} #3 \end{footnotesize}');
      fprintf(file.id, '%s\n', '  \begin{center} \texttt{File: #1} \end{center}');
      fprintf(file.id, '%s\n', '  \end{figure}');
      fprintf(file.id, '%s\n', '  \clearpage');
      fprintf(file.id, '%s\n', '}');
      fprintf(file.id, '%s\n', '\newcommand*{\inputSIDEFIG}[3]{');
      fprintf(file.id, '%s\n', '  \begin{sidewaysfigure}[t]');
      %       fprintf(file.id, '%s\n', '  \caption{#2}');
      fprintf(file.id, '%s\n', '  \centering\includegraphics[width=.75\textheight]{#1}');
      fprintf(file.id, '%s\n', '  \\ \noindent');
      fprintf(file.id, '%s\n', '  \begin{footnotesize} #3 \end{footnotesize}');
      % fprintf(file.id, '%s\n', '  \begin{center} \caption{#2} -- \texttt{File: #1} \end{center}');
      fprintf(file.id, '%s\n', '  \begin{center} \caption{#2} \end{center}');
      fprintf(file.id, '%s\n', '  \end{sidewaysfigure}');
      fprintf(file.id, '%s\n', '  \clearpage');
      fprintf(file.id, '%s\n', '}');

      fprintf(file.id, '\\lstset{numbers = left, escapechar = {}, numberstyle = \\tiny, basicstyle=\\ttlisting, breaklines = true, commentstyle = \\commentlisting, showstringspaces = true, keywordstyle = \\keywordlisting}\n'); %, char(163));
      fprintf(file.id, '%s\n', '\newfont{\ttlisting}{cmtt10 scaled 900}');
      fprintf(file.id, '%s\n', '\newfont{\commentlisting}{cmsl10 scaled 800}');
      fprintf(file.id, '%s\n', '\newfont{\keywordlisting}{cmb10 scaled 900}');

      fprintf(file.id, '%s\n', '\renewcommand*{\vec}{\boldsymbol}');

      fprintf(file.id, '%s\n', '\newcommand*{\inputTEXT}[2]{\section{#2}\lstinputlisting{#1}\vspace{2\baselineskip}}');
      fprintf(file.id, '%s\n', '\newcommand*{\inputDIARY}[2]{\part{#2}\lstinputlisting{#1}\vspace{2\baselineskip}}');
      fprintf(file.id, '%s\n', '\newcommand*{\inputLISTING}[2]{\section{#2}\lstinputlisting{#1}\vspace{2\baselineskip}}');
%       fprintf(file.id, '%s\n', '\newcommand*{\tttELMARmail}{\texttt{em@elmarmertens.ch}}');
%       fprintf(file.id, '%s\n', '\newcommand*{\linkELMARmail}{mailto:em@elmarmertens.ch}');
%       fprintf(file.id, '%s\n', '\newcommand*{\linkELMARmailRE}[1]{mailto:em@elmarmertens.ch?subject=#1}');
%       fprintf(file.id, '%s\n', '\newcommand*{\linkELMARmailDOC}{\linkELMARmailRE{Your\%20latexwrapper.m}}');
      fprintf(file.id, '\\date{%s}\n', datestr(now));

      %       if isempty(file.author)
      %          fprintf(file.id, '\\author{RESULTS STORED AT:\\thanks{Compiled with \\texttt{latexwrapper.m} from Elmar Mertens, email \\href{\\linkELMARmail}{\\tttELMARmail}. Download the program at \\href{http://www.elmarmertens.ch}{\\texttt{www.elmarmertens.ch}}.}\\\\ %s \\\\ %s}\n', ...
      %              strrep(unixfilename(file.dir), '$', '\$'), latexstr(file.name));
      %       else
      %          fprintf(file.id, '\\author{%s \\thanks{compiled with \\texttt{latexwrapper.m} from Elmar Mertens, email \\href{\\linkELMARmail}{\\tttELMARmail}. Download the program at \href{http://www.elmarmertens.ch}{\\texttt{www.elmarmertens.ch}.}}}\n', file.author);
      %       end


      if ~isempty(file.title)
         fprintf(file.id, '\\title{%s}\n', file.title);
      end

      % END OF PREAMBLE
      % LATEX DOC BODY STARTS HERE
      fprintf(file.id, '%s\n', '\begin{document}');
      if ~isempty(file.title)
         fprintf(file.id, '%s\n', '\maketitle');
      else
         fprintf(file.id, '%s\n', '%s', datestr(now));
      end

      fprintf(file.id, '%s\n', '\tableofcontents');
      fprintf(file.id, '%s\n', '\listoftables');
      fprintf(file.id, '%s\n', '\listoffigures');
      if ~isempty(file.pagestyle)
         fprintf(file.id, '\\pagestyle{%s}\n', file.pagestyle);
      end
      fprintf(file.id, '%s\n', '\newpage');

   case 'add'
      type2add = varargin{1};
      switch lower(type2add)
         case 'line'
            fprintf(file.id, '%s\n', varargin{2});
         case 'listing'
            file2add = varargin{2};
            fprintf(file.id, '\\lstset{language = MATLAB}\n');
            [jim, joe, jack] = fileparts(file2add);
            if any(strfind(jim, file.dir) == 1)
               joe = strcat(jim,'/',joe);
            end
            fprintf(file.id, '\\inputLISTING{%s}{%s}\n', file2add, strrep(strcat(joe,jack), '_', '\_'));
            fprintf(file.id, '\\lstset{language = {}}\n');
         case 'str2listing'
            str2add = varargin{2};
            fprintf(file.id, '\\lstset{language = MATLAB}\n');
            fprintf(file.id, '\\begin{lstlisting}\n');
            fprintf(file.id, '%s', str2add);
            fprintf(file.id, '\\end{lstlisting}\n');
         otherwise
            file2add = varargin{2};
            if length(varargin) > 2 && ~isempty(varargin{3})
               caption2add = varargin{3};
            else
               caption2add = file2add;
            end
            caption2add = regexprep(caption2add, '_', '\\_');
            if length(varargin) > 3
               comment2add = varargin{4};
            end
            switch lower(type2add)
               case {'tab', 'table'}
                  lcom = 'inputTAB';
               case 'sidetab'
                  lcom = 'inputSIDETAB';
               case {'fig', 'figure'}
                  lcom = 'inputFIG';
               case {'sidefig', 'sidewaysfigure'}
                  lcom = 'inputSIDEFIG';
               case {'diary'}
                  fprintf(file.id, '\\lstset{escapechar = {@}}\n');
                  lcom = 'inputDIARY';
               case {'text'}
                  lcom = 'inputTEXT';
               otherwise
                  error('type2add %s unknown', type2add)
            end
            if exist('comment2add', 'var')
               fprintf(file.id, '\\%s{%s}{%s}{%s}\n', lcom, file2add, caption2add, comment2add);
            else
               fprintf(file.id, '\\%s{%s}{%s}\n', lcom, file2add, caption2add);
            end
            fprintf(file.id, '\\newpage\n');
            if strcmpi(type2add, 'diary')
               fprintf(file.id, '\\lstset{escapechar = {}}\n');
            end
      end

   case 'dvi2pdf'
      if file.id ~= 0
         latexwrapper(file, 'compile');
      end
      if ispc % go via dvips because miktex's dvipdfm doesn't seem to support hyperref
         latexwrapper(file, 'makeps');
         latexwrapper(file, 'ps2pdf');
      else
         fprintf('dvipdfm %s.dvi ...', fullfile(file.dir, file.name));
         quietsystem(sprintf('cd %s && dvipdfm %s', file.dir, file.name))
         fprintf('done.\n');
      end
   case 'ps2pdf'
      if file.id ~= 0
         latexwrapper(file, 'compile');
      end
      if ~exist(fullfile(file.dir, strcat(file.name, '.ps')), 'file')
         latexwrapper(file, 'makeps');
      end
      lwd = cd;
      cd(file.dir); % cannot call shell with "cd file.dir && ..." in windows
      fprintf('ps2pdf %s.ps ...', fullfile(file.dir, file.name));
      quietsystem(sprintf('ps2pdf %s.ps', file.name))
      fprintf('done.\n');
      cd(lwd)
   case 'pdflatex'
      if file.id ~= 0
         file = latexwrapper(file, 'compilepdf');
      end
      wd = pwd;
      if isfield(file, 'dir')
         cd(file.dir)
      end
      fprintf('pdflatex %s.tex ...', fullfile(file.dir, file.name));
      quietsystem(sprintf('pdflatex -interaction=nonstopmode %s', file.name));
      fprintf('done.\n');
      cd(wd);
      
   case 'makeps'
      if file.id ~= 0
         latexwrapper(file, 'compile');
      end
      lwd = cd;
      cd(file.dir); % cannot call shell with "cd file.dir && ..." in windows
      fprintf('dvips %s.dvi ...', fullfile(file.dir, file.name));
      quietsystem(sprintf('dvips -q -z %s', file.name))
      fprintf('done.\n');
      cd(lwd)
   case {'compile'}

      if file.id ~= 0
         fprintf(file.id, '%s\n', '\end{document} ');
         file.id = fclose(file.id);
      end


      % switch diary off, unless it is the standard diary
      if strcmpi(get(0, 'Diary'), 'on') && ~strcmpi(get(0, 'DiaryFile'), fullfile(localroot, 'matlab.log'))
         diary off
      end

      % delete pdf and ps to avoid collisions
      wd = pwd;
      if isfield(file, 'dir')
         cd(file.dir)
      end
      if exist(strcat(file.name, '.pdf'), 'file')
         delete(strcat(file.name, '.pdf'))
      end
      if exist(strcat(file.name, '.ps'), 'file')
         delete(strcat(file.name, '.ps'))
      end
      [status, log] = system(sprintf('latex -interaction=nonstopmode %s', file.name));
      if status
         disp(log)
      else
         fprintf('LaTeX Compilation of\n\t "%s.tex"\nhas been successful (first run). Performing reruns ... ', fullfile(file.dir, file.name));
         quietsystem(sprintf('latex -interaction=nonstopmode %s', file.name));
         quietsystem(sprintf('latex -interaction=nonstopmode %s', file.name));
         fprintf('done.\n');
      end
      cd(wd);
      latexwrapper(file, 'clean')
   
    case {'compiledvi2pdf'}

      if file.id ~= 0
         fprintf(file.id, '%s\n', '\end{document} ');
         file.id = fclose(file.id);
      end


      % switch diary off, unless it is the standard diary
      if strcmpi(get(0, 'Diary'), 'on') && ~strcmpi(get(0, 'DiaryFile'), fullfile(localroot, 'matlab.log'))
         diary off
      end

      % delete pdf and ps to avoid collisions
      wd = pwd;
      if isfield(file, 'dir')
         cd(file.dir)
      end
      if exist(strcat(file.name, '.pdf'), 'file')
         delete(strcat(file.name, '.pdf'))
      end
      if exist(strcat(file.name, '.ps'), 'file')
         delete(strcat(file.name, '.ps'))
      end
      [status, log] = system(sprintf('latex -interaction=nonstopmode %s', file.name));
      if status
         disp(log)
      else
         fprintf('LaTeX Compilation of\n\t "%s.tex"\nhas been successful (first run). Performing reruns ... ', fullfile(file.dir, file.name));
         quietsystem(sprintf('latex -interaction=nonstopmode %s', file.name));
         quietsystem(sprintf('latex -interaction=nonstopmode %s', file.name));
         fprintf('done.\n');
      end
   
      fprintf('Converting DVI to PS  ... ');
      [status, log] = system(sprintf('dvips %s', file.name));
      if status
          disp(log)
      end
      fprintf('done.\n');
      
      fprintf('Converting PS to PDF  ... ');
      
      [status, log] = system(sprintf('ps2pdf %s.ps', file.name));
      if status
          disp(log)
      end
      fprintf('done.\n');

      cd(wd);
      latexwrapper(file, 'clean')

    case {'compilepdf'}

      if file.id ~= 0
         fprintf(file.id, '%s\n', '\end{document} ');
         file.id = fclose(file.id);
      end


      % switch diary off, unless it is the standard diary
      if strcmpi(get(0, 'Diary'), 'on') && ~strcmpi(get(0, 'DiaryFile'), fullfile(localroot, 'matlab.log'))
         diary off
      end

      % delete pdf to avoid collusions
      wd = pwd;
      if isfield(file, 'dir')
         cd(file.dir)
      end
%       if exist(strcat(file.name, '.pdf'), 'file')
%          delete(strcat(file.name, '.pdf'))
%       end
%       if exist(strcat(file.name, '.ps'), 'file')
%          delete(strcat(file.name, '.ps'))
%       end
      [status, log] = system(sprintf('pdflatex -interaction=nonstopmode %s', file.name));
      if status
         disp(log)
      else
         fprintf('PDFLaTeX Compilation of\n\t "%s.tex"\nhas been successful (first run). Performing reruns ... ', fullfile(file.dir, file.name));
         quietsystem(sprintf('pdflatex -interaction=nonstopmode %s', file.name));
%          quietsystem(sprintf('pdflatex -interaction=nonstopmode %s', file.name));
         fprintf('done.\n');
      end
      cd(wd);
      latexwrapper(file, 'clean')

   case 'edit'
      latexwrapper(file, 'open', 'tex');
   case 'log'
      latexwrapper(file, 'open', 'log');
   case 'open'

      if isempty(varargin)
         if ispc
            howopen = 'ps';
         else
            howopen = 'pdf';
         end
      else
         howopen = lower(varargin{1});
      end

      fprintf('Opening %s.%s ... \n', fullfile(file.dir, file.name), howopen)

      switch howopen
         case 'ps'
            if ~exist(fullfile(file.dir, strcat(file.name, '.ps')), 'file')
               latexwrapper(file, 'makeps')
            end
            if ispc
               dos(sprintf('start %s.%s', fullfile(file.dir, file.name), howopen));
            else
               unix(sprintf('kghostview %s.ps &', fullfile(file.dir, file.name)));
            end
         case 'pdf'
            if ~exist(fullfile(file.dir, strcat(file.name, '.pdf')), 'file')
               latexwrapper(file, 'compilepdf')
            end
            if ispc
               dos(sprintf('start %s.%s', fullfile(file.dir, file.name), howopen));
            else
               unix(sprintf('acroread %s.pdf &', fullfile(file.dir, file.name)));
            end
         case 'kpdf'
            if ~exist(fullfile(file.dir, strcat(file.name, '.pdf')), 'file')
               latexwrapper(file, 'compilepdf')
            end
            if ispc
               dos(sprintf('start %s.pdf', fullfile(file.dir, file.name)));
            else
               unix(sprintf('kpdf %s.pdf &', fullfile(file.dir, file.name)));
            end
         case 'acroread'
            if ~exist(fullfile(file.dir, strcat(file.name, '.pdf')), 'file')
               latexwrapper(file, 'compilepdf')
            end
            if ispc
               dos(sprintf('start %s.pdf', fullfile(file.dir, file.name)));
            else
               unix(sprintf('acroread %s.pdf &', fullfile(file.dir, file.name)));
            end
         case 'dvi'
            if ispc
               dos(sprintf('start %s.%s', fullfile(file.dir, file.name), howopen));
            else
               unix(sprintf('xdvi %s.dvi &', fullfile(file.dir, file.name)));
            end
         case 'xdvi'
            if ispc
               dos(sprintf('start %s.dvi', fullfile(file.dir, file.name)));
            else
               unix(sprintf('xdvi %s.dvi &', fullfile(file.dir, file.name)));
            end
         case 'tex'
            if ispc
               dos(sprintf('start %s.%s', fullfile(file.dir, file.name), howopen));
            else
               unix(sprintf('kwrite %s.pdf &', fullfile(file.dir, file.name)));
            end
         case {'emacs', 'xemacs'}
            system(sprintf('xemacs %s.tex &', fullfile(file.dir, file.name)));
         case 'log'
            type(fullfile(file.dir, strcat(file.name, '.log')))
         case {'shell', 'konsole'}
            if isunix
               unix(sprintf('konsole --workdir %s &', file.dir));
            else
               warning('em:msg', 'Opening with "%s" only supported on unix systems', howopen)
            end
         otherwise
            error('opening latex output with "%s" unknown.', howopen)
      end

   case 'clean'
      wd = pwd;
      if isfield(file, 'dir')
         cd(file.dir)
      end
      tmpfiles = {'*.tmp'};
      delete(tmpfiles{:})

      latexfiles = {'aux', 'toc', 'bbl', 'log', 'out', 'lof', 'lot', 'blg', 'tmp', 'touch'};
      for l = 1 : length(latexfiles)
         delfile = strcat(file.name, '.', latexfiles{l});
         if exist(delfile, 'file')
            delete(delfile)
         end
      end

      cd(wd);

   case 'script2compile' % create a script for compilation (useful if LaTeX cannot be called from Matlab)

      scriptname = fullfile(file.dir, strcat(file.name, '.sh'));
      [sid, message] = fopen(scriptname, 'wt');
      if sid < 0
         error('Could not open script file: "%s" message: "%s"', filename, message)
      end
      fprintf(sid, '%s\n', '');
      fprintf(sid, '#! /bin/bash');
      fprintf(sid, '# script to compile %s\n', file.name);
      fprintf(sid, 'latex %s ', file.name);
      fprintf(sid, '&& latex %s ', file.name);
      fprintf(sid, '&& latex %s ', file.name);
      fprintf(sid, '&& dvi2pdf %s ', file.name);
      fprintf(sid, '&& (acroread %s.pdf &) ', file.name);
      fprintf(sid, '&& echo && echo DONE\n');
      fprintf(sid, 'if [ -e %s.TOCOMPILE ] \nthen \nrm %s.TOCOMPILE \nfi\n', file.name, file.name);
      
      unix(sprintf('cd %s; touch %s.TOCOMPILE', file.dir, file.name));
      
   case {'close', 'stop'}

      if file.id
         fprintf(file.id, '%s\n', '\end{document} ');
         file.id = fclose(file.id);
      end
      wd = pwd;
      if isfield(file, 'dir')
         cd(file.dir)
      end

      % delete previous outputs to avoid collisions
      latexfiles = {'pdf', 'ps', 'dvi', 'touch'};
      for l = 1 : length(latexfiles)
         delfile = strcat(file.name, '.', latexfiles{l});
         if exist(delfile, 'file')
            delete(delfile)
         end
      end
      cd(wd)
   otherwise
      error('Command <<%s>> no recognized', command);
end

if nargout > 0
   outfile = file;
end

% -------------------------------------
function quietsystem(s)
% -------------------------------------

[status, log] = system(s);
if status
   disp(log)
end
