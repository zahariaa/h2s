function out = cellinsideout(in)
% cellinsideout: inverts the cell structure so, for example, if in is a 2x4
%    cell with each element being itself a 3x5 cell, out will be a 3x5 cell
%    with 2x4 cell elements.

szouter = size(in);
szinner = size(in{1});

out = repmat({cell(szouter)},szinner);
for si = 1:prod(szinner)
   for so = 1:prod(szouter)
      out{si}{so} = in{so}{si};
   end
end
return

