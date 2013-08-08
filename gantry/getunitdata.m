function ConvertStruct = getunitdata( mtype )
% GETUNITDATA internal function retrieving data structure for unit
% conversion. 


%error(nargoutchk(0,1, nargout,'struct'));

% load unit conversion data structure
ConvertStruct = aeroconvertdata(mtype);
