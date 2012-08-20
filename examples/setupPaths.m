%SETUPPATHS  Add paths for the Closest Point Method
%
%   TODO: need to make sure this works on Mac and Windows too.

disp('Adding paths for cp_matrices and surfaces');

% this is the logic of what we want to do but this is unix
% specific and doesn't work if you change directories
%addpath('../cp_matrices');
%addpath('../surfaces');
%addpath('../surfaces/tri');
%addpath('../surfaces/readply');
%addpath('../surfaces/tri2cp');

% TODO: probably only works if pwd is the the examples/ dir?
[cproot,name] = fileparts(pwd);

% add the cp_matrices directory
if isempty(strfind(path, 'cp_matrices'))
  addpath(fullfile(cproot, 'cp_matrices'));
end

% add the surfaces directory and subdirs
if isempty(strfind(path,'surfaces'))
  % this would add all subdirs
  %path(genpath(fullfile(cproot, 'surfaces')))
  addpath(fullfile(cproot, 'surfaces'));
  addpath(fullfile(cproot, 'surfaces', 'tri'));
  addpath(fullfile(cproot, 'surfaces', 'readply'));
  addpath(fullfile(cproot, 'surfaces', 'tri2cp'));
end



