%============================================================================
% Copyright (C) 2014, Heikki Hyyti
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%============================================================================

%set legend location to NorthWest, NorthEast, SouthWest or SouthEast where
%it least conflicts the data in plot. data is given in a matrix of
%different data sets in different columns. equally sampled data is assumed
%from first to last item from left to right

function location = legendLocation(data)
    dataLen = size(data,1);
    dataMaxVal = max(data(:));
    dataMinVal = min(data(:));
    first_part = data(1:ceil(dataLen/4),:);
    firstMaxVal = max(first_part(:));
    firstMinVal = min(first_part(:));
    last_part = data(floor(3*dataLen/4):end,:);
    lastMaxVal = max(last_part(:));
    lastMinVal = min(last_part(:));  
    
    selectLabels = {'NorthWest', 'NorthEast', 'SouthWest', 'SouthEast'};
    selectVals = zeros(4,1);
    selectVals(1) = dataMaxVal - firstMaxVal;
    selectVals(2) = dataMaxVal - lastMaxVal;
    selectVals(3) = firstMinVal - dataMinVal;
    selectVals(4) = lastMinVal - dataMinVal;
    
    [~,idx] = max(selectVals);
    location = selectLabels{idx};
end