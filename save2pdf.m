%SAVE2PDF Saves a figure as a properly cropped pdf
%
%   save2pdf(pdfname,directory,handle,dpi)
%
%   - pdfname: file name of pdf
%   - directory (optional): directory/destination directory of the output pdf file
%   - handle:  (optional) Handle of the figure to write to a pdf.  If
%              omitted, the current figure is used.  Note that handles
%              are typically the figure number.
%   - dpi: (optional) Integer value of dots per inch (DPI).  Sets
%          resolution of output pdf.  Note that 150 dpi is the Matlab
%          default and this function's default, but 600 dpi is typical for
%          production-quality.
%
%   Saves figure as a pdf with margins cropped to match the figure size.

%   (c) Gabe Hoffmann, gabe.hoffmann@gmail.com
%   Written 8/30/2007
%   Revised 9/22/2007
%   Revised 1/14/2007
%   Last modified by Hause Lin 19-10-20 21:21 hauselin@gmail.com

function save2pdf(pdfname,directory,handle,dpi)

if nargin<1 % if no pdfname provided, save as _unnamed.pdf in working directory
    pdfname = '_unnamed.pdf';
end
if ~contains(pdfname,'.pdf') % if file name doesn't contain .pdf, add it
    pdfname = [pdfname '.pdf'];
end
if nargin<2 % if no directory provided, save in current directory
    directory = '.';
end
if ~exist(directory)
    mkdir(directory)
end
pdfname = fullfile(directory,pdfname);
if nargin<3
    handle = gcf;
end
if nargin<4
    dpi = 300;
end

% Backup previous settings
prePaperType = get(handle,'PaperType');
prePaperUnits = get(handle,'PaperUnits');
preUnits = get(handle,'Units');
prePaperPosition = get(handle,'PaperPosition');
prePaperSize = get(handle,'PaperSize');

% Make changing paper type possible
set(handle,'PaperType','<custom>');

% Set units to all be the same
set(handle,'PaperUnits','inches');
set(handle,'Units','inches');

% Set the page size and position to match the figure's dimensions
paperPosition = get(handle,'PaperPosition');
position = get(handle,'Position');
set(handle,'PaperPosition',[0,0,position(3:4)]);
set(handle,'PaperSize',position(3:4));

% Save the pdf (this is the same method used by "saveas")
print(handle,'-dpdf',pdfname,sprintf('-r%d',dpi))

% Restore the previous settings
set(handle,'PaperType',prePaperType);
set(handle,'PaperUnits',prePaperUnits);
set(handle,'Units',preUnits);
set(handle,'PaperPosition',prePaperPosition);
set(handle,'PaperSize',prePaperSize);
