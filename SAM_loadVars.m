%% SAM_loadVars
% Load variables from a SAM NetCDF output file
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/17/2016
%
%%% Syntax
%   dat = SAM_getState(vars, fileloc)
%
%   dat = SAM_getState(vars, fileloc, castToDouble)
%
%%% Description
% Loads variables from a SAM NetCDF file. The function can return either a
% single variable in a matrix if only a single variable is requested or a
% struct containing matrices for each variable as fields if multiple
% variables are requested.
%
% Although this function was written for working with SAM output, it is
% general enough to work with most NetCDF files.
%
%
%%% Input Arguments
% *vars - variable name(s).*
% Either a single variable name or a cell array containing several variable
% names can be provided. The variable names must variable names in the
% NetCDF file. If a cell array with a single variable name is provided, the
% function will return a struct with a single field.
%
% *fileloc - path to the NetCDF file.*
% A single path must be provided. This function cannot read from multiple
% files in a single call.
%
% *castToDouble - convert numeric types to double precision (optional).*
% If this flag is set to 1, all numeric variables will be converted to
% double precision when they are read in.
%
%%% Output Arguments
% *dat - requested variables from the NetCDF file.*
% Either a matrix containing the data corresponding to a single variable
% (if a single string is provided as a variable name) or a structure
% containing variable data as fields (if a cell array of variable names is
% provided. The fields of the struct are named using the provided variable
% names converted to lowercase.
%
%%% <../test/html/SAM_loadVars_test.html Tests>
function dat = SAM_loadVars(vars, fileloc, varargin)

    % Parse input
    p = inputParser;
    checkVars = @(x) assert(iscell(x) || ischar(x),...
        'Variables must be given as a string or a cell array');
    checkFileloc = @(x) assert(ischar(x),...
        'File location must be a string');
    checkCast = @(x) assert(any(x == [0, 1]), ...
        'Double cast flag must be either 0 or 1');
    addRequired(p, 'vars', checkVars);
    addRequired(p, 'fileloc', checkFileloc);
    addOptional(p, 'castToDouble', 0, checkCast); 
    parse(p, vars, fileloc, varargin{:});
    castToDouble = p.Results.castToDouble;

    % Print information
%     fprintf(1,'Reading from NetCDF file: %s \n', fileloc);
    vn=0;
    vnmax=length(vars);
%     fprintf(1,'reading in variable #: %02d of %02d', vn, vnmax); 

    % Open the NetCDF file
    ncid3d = netcdf.open(fileloc, 'NC_NOWRITE');

    % Read in variables
    % Case 1: only a single variable
    if ~iscell(vars)
        vn=vn+1;
%         fprintf(1,'\b\b\b\b\b\b\b\b%02d of %02d', vn, vnmax); 
        ncname = char(vars);
        dat = netcdf.getVar(ncid3d, netcdf.inqVarID(ncid3d, ncname));
        if castToDouble && isnumeric(dat) 
            dat = double(dat);
        end
        
    % Case 2: multiple variables
    else
        dat=struct;
        for i=1:length(vars)
            vn=vn+1;
%             fprintf(1,'\b\b\b\b\b\b\b\b%02d of %02d', vn, vnmax); 
            vname = char(lower(vars(i)));
            ncname = char(vars(i));
            dat.(vname) = ...
                netcdf.getVar(ncid3d,netcdf.inqVarID(ncid3d,ncname));
            if castToDouble && isnumeric(dat.(vname))
                dat.(vname) = double(dat.(vname));
            end
        end
        fprintf(1,'\n');
    end

    % Close NetCDF file
    netcdf.close(ncid3d);

end