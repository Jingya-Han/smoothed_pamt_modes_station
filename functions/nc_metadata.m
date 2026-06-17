% Add global attributes (metadata)
function [ncfile1]=nc_metadata(ncfile1)
% 4 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% input--------------------------------------------------------------------
% ncfile1: name of ncfile
% output (does not need it actually)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncwriteatt(ncfile1, '/', 'author',      'Jingya Han');
ncwriteatt(ncfile1, '/', 'email',       'jh2423@cornell.edu');
ncwriteatt(ncfile1, '/', 'institution', 'Cornell University');
ncwriteatt(ncfile1, '/', 'reference',   'Han et al. (2026), Journal Name');
ncwriteatt(ncfile1, '/', 'history',     ['Created on ' datestr(now)]);

end