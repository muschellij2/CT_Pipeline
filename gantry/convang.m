function out  = convang(  in, uin, uout )
%  CONVANG Convert from angle units to desired angle units.
%   A = CONVANG( V, UI, UO ) computes the conversion factor from
%   specified input angle units, UI, to specified output angle units, UO, and
%   applies the conversion factor to the input, V, to produce the output,
%   A, in the desired units.  V and A are floating point arrays of size
%   M-by-N.	All of the values in V must have the same unit conversions from
%   UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'deg'         :degrees       
%      'rad'         :radians
%      'rev'         :revolutions
%
%   Example:
%
%   Convert three angles from degrees to radians:
%      a = convang([3 10 20], 'deg','rad')
%
%   See also CONVACC, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVFORCE, CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP, CONVVEL.

%   Copyright 2000-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/11/01 23:38:52 $

if ~isfloat( in )
    error('aero:convang:notfloat','Angle vector is not a floating point number');
end

slope = unitconversion('angle conversion',uin,uout,0);

out = in.*slope;
