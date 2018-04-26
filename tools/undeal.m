function varargout = undeal(outputnums,fcn)
% undeal.m: Wrapper to select output(s) for a given anonymous function
% 
% Ex. equivalent calls
% [center,centerCI,radius,radiusCI] =estimateHypersphere(X);
% [center,radius] = undeal([1 3],@() estimateHypersphere(X));
% cntr_radius_cell= undeal({1 3},@() estimateHypersphere(X));
% 
% Ex. cross-validate 3rd argument
% cv = cvindex(size(X,1),10);
% radiusEsts = cv.crossvalidate(@(X) undeal(3,@() estimateHypersphere(X)),X);
%
% 2018-04-26 AZ Created

%% PRELIMINARIES!
if iscell(outputnums),   CELLOUT = true;   outputnums = [outputnums{:}];
else                     CELLOUT = false;
end
maxo = max(outputnums);

%% EVALUATE!
eval(sprintf('[%svarargout{%u}] = fcn();',sprintf('varargout{%u},',1:maxo-1),maxo));

%% DEAL!
if CELLOUT,  varargout{1} = varargout(outputnums);
else         varargout    = varargout(outputnums);
end

return

