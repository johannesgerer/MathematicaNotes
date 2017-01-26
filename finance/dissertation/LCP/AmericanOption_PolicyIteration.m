% (c) Jan Hendrik Witte, Mathematical Institute, University of Oxford
% witte@maths.ox.ac.uk

% Matlab Code for pricing an American Put in a Standard Black-Scholes Model
% Eg see "The Mathematics of Financial Derivatives" by Wilmott, Howison, Dewynne)
% for the general model
% Policy Iteration as in Algorithm 2.1 in http://arxiv.org/abs/1012.4976

function [S,V] = AmericanOption_PolicyIteration(M, N, K, r, sigma, T, theta, DoPlot)
% eg run: AmericanOption_PolicyIteration(100, 100, 1, .03, .6, 1, .5, 1)

% M number of time steps
% N number of space steps 
% K strike
% T time horizon
% theta - time weighting parameter
% theta = 0, 1/2, 1 is fully explicit, crank-nicolson, fully implicit
% DoPlot - set 1 to show plot, else 0

%
% setup time and space steps
%
Smax = 6 * K;
h    = Smax / N;
k    = T / M;
n    = [0:N]';
S    = n*h;

%
% setup coefficient matrix (implicit part)
%
a = -0.5*theta*k*(sigma^2*n.^2 - r*n);
b = 1 + theta*k*(sigma^2*n.^2 + r);
c = -0.5*theta*k*(sigma^2*n.^2 + r*n);

mat = spdiags([c b a],[-1:1],N+1,N+1)';

%
% setup coefficient matrix (explicit part)
%
A = 0.5*(1 - theta)*k*(sigma^2*n.^2 - r*n);
B = 1 - (1 - theta)*k*(sigma^2*n.^2 + r);
C = 0.5*(1 - theta)*k*(sigma^2*n.^2 + r*n);

MAT = spdiags([C B A],[-1:1],N+1,N+1)';

% intial condition
payoff = max(K - S, 0);
NormPf = norm(payoff,inf);
V      = payoff;

% parameters for american solver
tol          = 1e-6;
kmax         = 25;
E            = speye(N+1);

% time stepping
for m = [0:M-1]
    RHS     = MAT * V;
    normRHS = norm(RHS,inf);
    Vk      = V;
    % start policy iteration
    err1    = -2*tol; err2 = 2*tol;
    k       = 0; 
    while (((err1 < -tol) || (err2 < -tol) || (err3 > tol))) && (k<kmax)
        Idx            = (mat*Vk-RHS) > (Vk-payoff);
        mattemp        = mat;
        RHStemp        = RHS;
        mattemp(Idx,:) = E(Idx,:);
        RHStemp(Idx)   = payoff(Idx);
        % solve linear system
        Vnew = mattemp\RHStemp;
        % check accuracy
        LHS  = mat*Vnew;
        err1 = min(LHS - RHS)/normRHS;
        err2 = min(Vnew - payoff)/NormPf;
        err3 = max(LHS(Vnew>payoff + tol*NormPf) - RHS(Vnew>payoff + tol*NormPf))/normRHS;
        % prepare next iteration
        Vk   = Vnew;
        k    = k+1; 
    end 
    if (k>=kmax) disp(k); warning('kmax has been reached'); end
    V = Vk;  
end

% plot the t=0 value
if DoPlot ==1
    figure(1)
    plot(S,V,'LineWidth',1)
    title('American Put Option')
    xlabel('Asset Price')
    ylabel('Option Value')
end

