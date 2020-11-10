function [fval,grad] = maxRadiusGivenCenter(m,X)
% maxRadiusGivenCenter: objective function to optimize, with gradient

grad = X - repmat(m(:)',[size(X,1) 1]);
fval = sqrt(sum(grad.^2,2));

% Apply max function
[fval,ix] = max(fval);
grad = -grad(ix,:)/fval;

return
end
