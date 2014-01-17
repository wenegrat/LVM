function uTest_BlockMean(doSpeed)
% Automatic test: BlockMean M and Mex
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_BlockMean(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R0n V:026 Sum:bdBKcYc4PlY0 Date:22-Sep-2010 02:30:25 $
% $License: BSD $
% $File: Tools\UnitTests_\uTest_BlockMean.m $
% History:
% 022: 11-Mar-2010 22:30, Test rectangular blocks.
%      M-version of isEqualTol_l included! The test function published in the
%      FEX (21-Jul-2009) failed, because I forgot to add this function. It seems
%      like nobody run the test, or nobody cares for a failing test, or nobody
%      uses the function at all after downloading. Now this test function should
%      run and I'm aware that I've overestimated its importance for others.  :-|

% Initialize: ==================================================================
LF = char(10);

if nargin == 0
   doSpeed = true;
end
minT = eps;

% Do the work: =================================================================
% Hello:
disp(['==== Test BlockMean:  ', datestr(now, 0), LF, ...
      '  Version: ', which('BlockMean'), LF]);

TypeList = {'uint8', 'double'};
for iType = 1:length(TypeList)
   aClass = TypeList{iType};
   disp(['== Test input type [', upper(aClass), ']']);
   
   X = mycast(zeros(10, 10), aClass);
   failed = 0;
   try
      Y = BlockMean(X);
      failed = 1;
   catch
      Y = [];
      disp(['  ok: Bad input recognized:', LF, '    ', lasterr]);
   end
   if length(Y) || failed
      error(['*** ', mfilename, ': Too friendly for 1 input!']);
   end
   
   Y = BlockMean(mycast([], aClass), 3);
   if length(Y) || (isa(Y, aClass) == 0)
      error(['*** ', mfilename, ': Failed on empty matrix']);
   end
   disp('  ok: Empty matrix');
      
   for dim3 = 1:4
      if dim3 > 1
         dim3S = sprintf(' x %d', dim3);
      else
         dim3S = '';
      end
      
      X = mycast(ones(3, 3, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(1, 1, dim3))
         disp(['  ok: [3 x 3', dim3S, '] => [1 x 1', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [3 x 3', dim3S, ']']);
      end
      
      X = mycast(ones(4, 3, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(1, 1, dim3))
         disp(['  ok: [4 x 3', dim3S, '] => [1 x 1', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [4 x 3', dim3S, ']']);
      end
      
      X = mycast(ones(5, 3, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(1, 1, dim3))
         disp(['  ok: [5 x 3', dim3S, '] => [1 x 1', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [5 x 3', dim3S, ']']);
      end
      
      X = mycast(ones(6, 3, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(2, 1, dim3))
         disp(['  ok: [6 x 3', dim3S, '] => [2 x 1', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [6 x 3', dim3S, ']']);
      end
      
      X = mycast(ones(3, 4, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(1, 1, dim3))
         disp(['  ok: [3 x 4', dim3S, '] => [1 x 1', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [3 x 4', dim3S, ']']);
      end
      
      X = mycast(ones(3, 5, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(1, 1, dim3))
         disp(['  ok: [3 x 5', dim3S, '] => [1 x 1', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [3 x 5', dim3S, ']']);
      end
      
      X = mycast(ones(3, 6, dim3), aClass);
      Y = BlockMean(X, 3);
      if isequal(Y, ones(1, 2, dim3)) && strcmp(class(X), class(Y))
         disp(['  ok: [3 x 6', dim3S, '] => [1 x 2', dim3S, ']']);
      else
         error(['*** ', mfilename, ': Failed on [3 x 6', dim3S, ']']);
      end
      
      X = repmat(mycast(reshape(1:35, 5, 7), aClass), [1, 1, dim3]);
      for iX = 1:6
         for iY = 1:8
            Y = BlockMean(X, iX, iY);
            Z = localBlockMean(X, iX, iY);
            if ~isequal(Y, Z)
               error(['*** ', mfilename, ': BlockMean([5 x 7', dim3S, ...
                     '], ', sprintf('%d, %d', iX, iY), ')']);
            end
         end
      end
      disp(['  ok: (5 x 7', dim3S, '), 1..6, 1..8)']);
   end
   
   % ---------------------------------------------------------------------------
   fprintf(['\n==  Random tests (5 sec, ', upper(aClass), ') ...']);
   % Need surprisingly high limit for 64 bit precision:
   % 16 * 16 eps relative error for numbers up to 255
   Tol   = 2048 * eps;
   iTime = cputime;
   nLoop = 0;
   while cputime - iTime < 5
      nLoop = nLoop + 1;
      if mod(nLoop, 100) == 0
         fprintf('.');
      end
      
      dim = [fix(rand(1, 2) * 100), fix(rand(1, 4) * 3 + 1)];
      dim = dim(1:max(2, max(find(dim > 1))));  %#ok<MXFND> Matlab 6 compatible
      X   = mycast(fix(rand(dim) * 255), aClass);
      V   = fix(rand * 16) + 1;
      W   = fix(rand * 16) + 1;
      
      try
         Y1 = BlockMean(X, V, W);
      catch
         fprintf('\n');
         error(['*** ', mfilename, ': Test crashed:', char(10), lasterr]);
      end
      Y2 = localBlockMean(X, V, W);
      switch aClass
         case 'double'
            eq = isEqualTol_l(Y1, Y2, Tol);
         case 'uint8'
            eq = isequal(Y1, Y2);
      end
      if eq == 0
         error(['*** ', mfilename, ': Failed on test data.']);
      end
   end
   fprintf(' %d tests ok\n', nLoop);
   
   % ===========================================================================
   disp([char(10), '== Speed tests (', upper(aClass), ') ...']);
   
   % Find a suiting number of loops:
   X = mycast(fix(rand(100, 32, 3) * 255), aClass);
   if doSpeed
      iLoop     = 0;
      startTime = cputime;
      while cputime - startTime < 1.0
         a = [];   %#ok<*NASGU> % Important to release the memory
         a = BlockMean(X, 10);
         b = [];
         b = localBlockMean(X, 10);
         iLoop = iLoop + 1;
      end
      nLoops = 100 * ceil(iLoop / ((cputime - startTime) * 100));
      disp([sprintf('  %d', nLoops), ' loops on this machine.']);
   else
      disp('  Use at least 2 loops (displayed times are crazy!)');
      nLoops = 2;
   end
   
   % -----------------
   disp([char(10), '  100x32x3 array in chunks of 10x10 frames:']);
   tic;
   for i = 1:nLoops
      a = [];
      a = BlockMean(X, 10, 10);
   end
   cTime = toc + minT;
   drawnow;             % Allow checks for Ctrl-C
   
   tic;
   for i = 1:nLoops
      a = [];
      a = localBlockMean(X, 10, 10);
   end
   mTime = toc + minT;
   
   disp(['    local M-version: ', sprintf('%.2f', mTime), LF, ...
         '    BlockMean:       ', sprintf('%.2f', cTime), '  ==> ', ...
         sprintf('%.1f', 100 * cTime / mTime), '% of Matlab time', LF]);
   
   % -----------------
   disp('  100x32x3 in chunks of 5x5 frames:');
   tic;
   for i = 1:nLoops
      a = [];
      a = BlockMean(X, 5, 5);
   end
   cTime = toc + minT;
   drawnow;
   
   tic;
   for i = 1:nLoops
      a = [];
      a = localBlockMean(X, 5, 5);
   end
   mTime = toc + minT;
   
   drawnow;
   disp(['    local M-version: ', sprintf('%.2f', mTime), LF, ...
         '    BlockMean:       ', sprintf('%.2f', cTime), '  ==> ', ...
         sprintf('%.1f', 100 * cTime / mTime), '% of Matlab time', LF]);
   
   % -----------------
   mLoops = round(nLoops / 100);
   fprintf('  1024x768x3 in chunks of 4x4 frames:  (%d loops)\n', mLoops);
   X = mycast(rand(1024, 768, 3) * 255, aClass);
   
   tic;
   for i = 1:mLoops
      a = [];
      a = BlockMean(X, 4, 4);
   end
   cTime = toc + minT;
   drawnow;
   
   tic;
   for i = 1:mLoops
      a = [];
      a = localBlockMean(X, 4, 4);
   end
   mTime = toc + minT;
   
   drawnow;
   disp(['    local M-version: ', sprintf('%.2f', mTime), LF, ...
         '    BlockMean:       ', sprintf('%.2f', cTime), '  ==> ', ...
         sprintf('%.1f', 100 * cTime / mTime), '% of Matlab time', LF]);
   
end  % for TypeList

% Goodbye:
disp('BlockMean seems to work well.');

return;

% ******************************************************************************
function X = mycast(X, aClass)
switch aClass
   case 'double'
   case 'uint8'
      X = uint8(X);
   otherwise  % Programming error:
      error('*** TextBlockMean: Unknown type.');
end
return;

% ******************************************************************************
function Y = localBlockMean(X, V, W)
% Local version of BlockMean using SUM(SUM(X))
if nargin < 3
   W = V;
end
S = size(X);
M = S(1) - mod(S(1), V);
N = S(2) - mod(S(2), W);
if M * N == 0
   Y = X([]);  % Copy type of X
   return;
end
MV = M / V;
NW = N / W;

% Cut and reshape input such that the 1st and 3rd dimension have the lengths V
% and W:
XM = reshape(X(1:M, 1:N, :), V, MV, W, NW, []);

% Different methods depending on the type of the input:
if isa(X, 'double')
   Y = sum(sum(XM, 1), 3) .* (1.0 / (V * W));
elseif uint8(0.8) == 1  % UINT8 of Matlab7 rounds:
   % NOT .*1/(V*W) !!! Otherwise the rounding differs from C-Mex
   Y = uint8(sum(sum(XM, 1), 3) ./ (V * W));
else                    % UINT8 of Matlab6 truncates, so round manually:
   Y = uint8(round(sum(sum(XM, 1), 3) ./ (V * W)));
end

% Remove singleton dimensions:
S(1) = MV;
S(2) = NW;
Y    = reshape(Y, S);

return;

% ******************************************************************************
function Equal = isEqualTol_l(x, y, Tol)
% Compare two double arrays with absolute tolerance
% Equal = isEqualTol_l(x, y, Tol)
% If the maximal difference between two double arrays is smaller than Tol, they
% are accepted as equal.
% As in Matlab's ISEQUAL, NaNs are treated as inequal, so isEqualTol_l(NaN, NaN)
% is FALSE. If you need comparable NaNs use ISEQUALWITHEQUALNANS, although this
% name is horrible and it is not existing in Matlab5.3.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2007-2010 matlab.THISYEAR(a)nMINUSsimon.de

% Was JRev: R0c V:022 Sum:mP3jthZjqvhS Date:12-Mar-2010 01:41:33
% Was File: Tools\GLMath\isEqualTol_l.m

% Simplified version!

% Try the exact and fast comparison at first:
Equal = false;
if isequal(size(x), size(y))
   % Same as "if all(abs(xMy(:)) <= Tol)", but faster:
   xMy = x - y;
   if all(or((abs(xMy) <= Tol), (x == y)))   % is FALSE for NaNs
       Equal = true;
   end
end

return;
