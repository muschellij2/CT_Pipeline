function [ varargout ] = unitconversion( mtype, iunit, ounit, useport )
% UNITCONVERSION internal function containing common algorithm for unit
% conversion. 

error(nargoutchk(0,3, nargout,'struct'));

if ~ischar( mtype )
    error('aero:unitconversion:notchar1','Unit conversion type is not a string');
end

if ~ischar( iunit )
    error('aero:unitconversion:notchar2','First unit declaration is not a string');
end

if ~ischar( ounit )
    error('aero:unitconversion:notchar3','Second unit declaration is not a string');
end

if ~isnumeric( useport )
    error('aero:convlength:notnumber4','Port flag is not a number');
end

% load unit conversion data structure
ConvertStruct = getunitdata(mtype);

units = {ConvertStruct.tdata.unit};

% find index for input unit conversion data
inidx = find(strcmpi(iunit,units),1);

if isempty( inidx )
    error('aero:unitconversion:unknownunit','Unknown %s input unit, %s', mtype, iunit);
end

% find index for output unit conversion data
outidx = find(strcmpi(ounit,units),1);

if isempty( outidx )
    error('aero:unitconversion:unknownunit','Unknown %s output unit, %s', mtype, ounit);
end

slope_in = ConvertStruct.tdata(inidx).slope;
slope_out = ConvertStruct.tdata(outidx).slope;
varargout{1} = slope_in/slope_out;

if strcmpi(mtype,'temperature conversion')
    bias_in = ConvertStruct.tdata(inidx).bias;
    bias_out = ConvertStruct.tdata(outidx).bias;
    varargout{2} = ( bias_in - bias_out )/slope_out;
end

if useport
    ports(1).txt = ConvertStruct.tdata(inidx).unit;
    ports(2).txt = ConvertStruct.tdata(outidx).unit;
    if strcmpi(ConvertStruct.tdata(inidx).unit,'naut mi')
        ports(1).txt = 'n.mi';
    end
    if strcmpi(ConvertStruct.tdata(outidx).unit,'naut mi')
        ports(2).txt = 'n.mi';
    end
    if strcmpi(mtype,'temperature conversion')
        varargout{3} = ports;
    else
        varargout{2} = ports;
    end
end
