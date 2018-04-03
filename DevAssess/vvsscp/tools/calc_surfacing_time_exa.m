clear all
close all
clc

% % % required user input 
folder = '../r1/';              % path to mfix exa ascii files
fnamebase = 'R_DES_DATA_';      % base file name 
dt  = 5.0d-4;                   % time step    
int = 40;                       % ascii write interval 
T   = 500*int;                  % int of last file
Ni  = 4;                        % number of intruder particles 
Nc  = 4300;                     % (initial) number of common particles
% % % required user input 


tsurf = [-1 -1 -1 -1];    % init to -ve; only check if intruders have "surfaced" if -ve to avoid intruderes "diving" again after they've surfaced

for tt = 0:int:T
    % import data
    fname = strcat(folder,fnamebase,num2str(tt,'%05i'));
    fid = fopen(fname,'r');
    Np = fscanf(fid, '%i', 1);
    jnk = fscanf(fid, '%i', 4);
    pdat = zeros(Np,21);
    for ii = 1:Np
        pdat(ii,:) = fscanf(fid, '%f', 21);
    end
    fclose(fid);
    
    % calc y-cdf of common particles and check if intruders have surfaced
    yi = pdat(1:Ni,1);        % y pos of (i) #NOTE: in exa put long dim first
    yc = pdat(Ni+1:Np,1);     % y pos of (c) #NOTE: in exa put long dim first
    yc_sort = sort(yc);       % sort
    yc95 = yc_sort(4085);     % pick out the 95th-percentile (c) 
    for ii = 1:Ni 
        if tsurf(ii) < 0 
            if yi(ii) > yc95
                tsurf(ii) = double(tt)*dt;
            end
        end
    end
end

% sort tsurf and print
tmp = sort(tsurf); 
for ii = 1:Ni
    fprintf('%10.2e\n',tmp(ii)); 
end

