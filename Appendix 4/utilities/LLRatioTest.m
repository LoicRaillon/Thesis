function p = LLRatioTest(LL_large,LL_small,dof)

small   = abs(LL_small);
large   = abs(LL_large);
ratio   = -2*(small-large);
% p       = 1 - gammainc(ratio/2, dof/2);
p       = 1 - pchisq(ratio, dof);


