%ETOPO1LOAD
% Downloads and extracts the ETOPO1 data from NOAA website. Uses curl on
% UNIX and built-in websave otherwise.
%   See also IN_DOMAIN
s_old = pwd;
s = char(mfilename('fullpath'));
cd(s(1:end-10))                         % cd to the mfile folder; 
                                        % The constant 12 depends on the
                                        % length of the filename.
url = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/cell_registered/binary/etopo1_bed_c_i2.zip';
filename = '../Output/etopo1_bed_c_i2.zip';
fprintf('Please wait while the ETOPO data is downloaded from\n');
fprintf('https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO:1/data/bedrock\n');
if isunix
    !curl  -o ../Output/etopo1_bed_c_i2.zip 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/cell_registered/binary/etopo1_bed_c_i2.zip'
else
    websave(filename,url)
end
unzip('../Output/etopo1_bed_c_i2.zip','../Output/')
fprintf('\n\n... done.\n');
clc;
cd(s_old)