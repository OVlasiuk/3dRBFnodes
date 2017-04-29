%ETOPO1LOAD
% Downloads and extracts the ETOPO1 data from NOAA website.
%   See also IN_DOMAIN
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-10))                         % cd to the mfile folder; 
                                        % The constant 12 depends on the
                                        % length of the filename.

fprintf('Please wait while the ETOPO data is downloaded.\n');
urlwrite('https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/cell_registered/binary/etopo1_bed_c_i2.zip','../Output/etopo1_bed_c_i2.zip')
unzip('../Output/etopo1_bed_c_i2.zip','../Output/')
fprintf('\n... done.\n');
cd(s_old)