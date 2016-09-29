function [xk] = sart2(A,b,K,x0,options)
% Update on Apr.10th 2015
%   1. Also input A' to speed up the recon
% Gongting Wu, UNC Chapel Hill

% Check that at least 3 inputs are given.
if nargin < 3
    error('Too few input arguments')
end
% 
% const = const^2;
% A = sqrt(Ax + Ay*const);
% clear Ax Ay const

[m,n] = size(A);

% Compute the diag matrix
Apj = full(sum(abs(A),1));
Aip = full(sum(abs(A),2));
W = 1./Aip;
I = (W == Inf);
W(I) = 0;
V = 1./Apj';
I = (V == Inf);
V(I) = 0;
clear Apj Aip I

AT = A';

% Check that the sizes of A and b match.
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
end

% Default value for x0.
if nargin < 4 || isempty(x0)
    x0 = zeros(n,1);
end

% Check the size og x0.
if size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem')
end

% Default values.
if nargin < 5
    % There must be a maximum number of iterations.
    if isempty(K)
        error('No stopping rule specified')
    end 

    % Default is no nonnegativity or box constraint.
    nonneg = false;
    boxcon = false;

else
% Check the contents of options if present.
    
    % Nonnegativity.
    if isfield(options,'nonneg')
        nonneg = options.nonneg;
    else
        nonneg = false;
    end
    
    % Box constraints [0,L].
    if isfield(options,'box')
        nonneg = true;
        boxcon = true;
        L = options.box;
    else
        boxcon = false;
    end
    
       
end % end if nargin includes options.

% Initialize the values.
rxk = b - A*x0;
xk = x0;

for k=1:K
    
    xk = xk + options.lambda(k)*(V.*(AT*(W.*rxk)));
    
    % Nonnegativity and box constraints.
    if nonneg, xk(xk<0) = 0; end
    if boxcon, xk(xk>L) = L; end
    
    % New residual.
    rxk = b - A*xk;
    
end