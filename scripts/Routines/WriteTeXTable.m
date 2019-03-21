function [] = WriteTeXTable(fid, header, style, tableBody, aboveTabular, belowTabular)
%% WriteTeXTable - Writing tex tables from Matlab.
%
% NOTE: below, I will use NCOL and NROW to denote the maximal number of
% columns in the TeX table ('maximal' in the sense that I ignore
% multicolumn), while NROW is the number of rows of data, excluding
% header information.
%
% Input arguments
% ---------------
% FID         File ID you got before running this program by executing
%             something like
%
%               fid = fopen('TableName.tex', 'w')
%
%             Can also be NaN or something empty ([], '', {}) to print
%             to screen for testing purposes.
%
% HEADER      Cell array of column titles/headers (multiple rows
%             accepted)
%
%             NOTE: A cell array of column titles/headers can look like
%             this
%
%                  {'Col1', 'Col2', 'Col3', 'Col4'};
%             OR
%                  {'Multi1', NaN, NaN, 'Col4'};
%
%             where having a NaN after text will cause that text to
%             span its column and all subsequent NaN columns via
%             multicolumn. Therefore, only text entries show up as
%             column titles/headers, while NaNs can be used to signify
%             multi-column spanning. A cline will automatically be
%             placed under.
%
% STYLE       Latex style for the columns like 'r|cccc' or 'l|rr|rr|'.
%
% TABLEBODY   Array or cell matrix of table content to write. Though
%             headers are treated differently from data to display, no
%             such distinction is made for the first column or couple of
%             columns for labels. Just throw them in here.
%
%             TABLEBODY as numeric araray
%             ---------------------------
%             The numbers will be displayed to significant digits
%             determined by this program. If you would like to display
%             them under a different choice of significant digits,
%             store the values as strings in a cell array and see the
%             next section.
%
%             TABLEBODY as cell array
%             -----------------------
%             This gives the most control, and different entries of
%             TABLEBODY will be treated differently depending upon the
%             data type.
%
%             In particular, here is the treatment by data type:
%             - TEXT      Write as is
%             - A NUMBER  Convert to significant digits determined by this
%                         program
%             - NaN       Take the previous value in the row, and
%                         stretch it via multicolumn into this row too.
%                         Multiple NaNs will stretch across more columns
%             - EMPTY     Write nothing
%
% ABOVETABULAR  What to write in between \begin{table} and
%               \begin{tabular}{style}. Can be a string or a cell of
%               strings (which would be written one entry per line).
% BELOWTABULAR  What to write in between \end{tabular} and \end{table}.
%               Can be a string or a cell of strings (which would be
%               written one entry per line).
%
% SIMPLE EXAMPLE
% --------------
%
%   row_names  = {'row1'; 'row2'; 'row3'};
%   data2write = randn(3);
%   header     = {'col1', 'col2', 'col3'};
%   table_data = [row_names, num2cell(data2write)];
%   style      = 'r|ccc';
%
%   fid = fopen('file.tex', 'w');
%   write_table_tex(fid, header, style, table_data, 'Random Numbers')
%   fclose(fid);
%

%% Define some helper functions

  % Inline if statement
  iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

  % Function to format numbers with specific significant digits and to
  % format text into TeX friendly mode (& -> \&)
  texify = @(x) x;
  fmt = @(x) iif(isnumeric(x) && ~isempty(x) && abs(x) < 1,  sprintf('%.3f', x), ...
                 isnumeric(x) && ~isempty(x) && abs(x) < 10, sprintf('%.2f', x), ...
                 isnumeric(x),                               sprintf('%.1f', x), ...
                 isempty(x),                                 '', ...
                 true,                                       @() texify(x));

  % To write a single line to the file, writeline(string) rather than
  % fprintf(fid, string). Also, no need to escape special chars this way
  if isempty(fid) || isnan(fid)
    writeline = @(s) fprintf('%s\n', s);
  else
    writeline = @(s) fprintf(fid, '%s\n', s);
  end
  writerow  = @(rowcell) writeline([strjoin(rowcell, ' & '), '\\']);

  % Write each element of a cell as a separate line;
  writeCell = @(towrite) cellfun(@(txt) writeline(fmt(txt)), towrite);

  % Takes either a string and writes that to a line, or a cell and
  % writes each entry as a separate line
  writeStringCell = @(towrite) iif( ischar(towrite), @() writeCell({towrite}), ...
                                    iscell(towrite), @() writeCell( towrite ));


%% Write the beginning of the table

  writeline('\begin{table}[htpb!]');
  if exist('aboveTabular', 'var') && ~isempty(aboveTabular)
    writeStringCell(aboveTabular);
  end
  writeline('\centering');
  writeline(['\begin{tabular}{', style, '}']);

% Get the columns after which we have a | vertical line (with 0 being a
% possibility if there is a vertical line before col 1 (so after row 0)

  vertLines = [];
  col = 0;
  for n = 1:length(style)
    if style(n) == '|'
      vertLines(end+1) = col;
    else
      col = col + 1;
    end
  end

%% Write header, then table body

  toWrite          = {header, tableBody};
  useClines        = [1 0];
  doubleHlineAfter = [1 0];

  % Loop over objects to write
  for w = 1:length(toWrite)
    writing = toWrite{w};
    [Nrow, Ncol] = size(writing);

    % Loop over rows within that object
    for r = 1:Nrow
      row = writing(r,:);
      if ~iscell(row)
        row = num2cell(row);
      end

      % Loop over columns of the row and mark its non-NaN entries
      notNaN = zeros(1, Ncol);
      for n = 1:Ncol
        if isempty(row{n}) || ischar(row{n}) || ~isnan(row{n})
          notNaN(n) = 1;
        end
      end
      notNaN = find(notNaN);

      % Loop over columns in row with actual text (rather than NaN, parse
      % text, account for multicolumns, and store in cell for row to be
      % written
      rowEntries = cell(1,length(notNaN));
      clines     = '';
      for n = 1:length(notNaN)

        % Start at current entry
        startCol = notNaN(n);

        % Stop 1 before next entry or at end of row
        if n < length(notNaN)
          stopCol = notNaN(n+1)-1;
        else
          stopCol = Ncol;
        end

        % Store multicolumn or single column
        colGap = stopCol- startCol;
        if colGap >= 1
          % Check for left and right verts
          leftVert = ''; rightVert = '';
          if any(vertLines == (startCol-1)), leftVert  = '|'; end
          if any(vertLines == stopCol),      rightVert = '|'; end

          rowEntries{n} = sprintf('\\multicolumn{%d}{%sc%s}{%s}', ...
                                  colGap+1, leftVert, rightVert, fmt(row{startCol}));
          clines = [clines sprintf('\\cline{%d-%d}', startCol, stopCol)];
        else
          rowEntries{n} = fmt(row{startCol});
        end
      end  % End row construction

      % Write row and clines
      writerow(rowEntries);
      if useClines(w) && ~isempty(clines)
        writeline(clines);
      end
    end

    if doubleHlineAfter(w)
      writeline('\hline\hline');
    end

  end % End header and body writing

% Conclude
writeline('\hline');
writeline('\end{tabular}');

% Caption
if exist('belowTabular', 'var') && ~isempty(belowTabular)
  writeStringCell(belowTabular);
end
writeline('\end{table}');

writeline('');

end
