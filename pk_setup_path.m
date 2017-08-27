subdirs = {'.', 'util', 'ops', 'kernels'};
nSubdirs = size(subdirs, 2);
for k = 1 : nSubdirs
  p = sprintf('%s/%s', pwd, subdirs{k});
  addpath(p);
end
clear p k nSubdirs subdirs;
global pk_home_dir;
pk_home_dir = pwd;
