function tf = isoctave()
%ISOCTAVE   Returns true if running GNU Octave
%   Example:
%   if isoctave()
%     octave_specific_code...
%   else
%     matlab_specific_code...
%   end
%
%   There seems to be some small overhead for this so probably not
%   a good idea to use it in a tight loop.

  tf = exist('octave_config_info', 'builtin');

