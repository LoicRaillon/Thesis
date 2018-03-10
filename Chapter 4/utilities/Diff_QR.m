% QR rectangular matrices differentiation 
% 
% Modern square-root Kalman filter implementations use QR factorization 
% at each iteration step of the filter: given a matrix A (pre-array), 
% assume that we apply an orthogonal operator Q to the A (i.e. QA=R) 
% so as to get the upper-triangular post-array R.
% In general, the A and R has the following structure: 
% one iteration of the array SR filters has the form
%    Q [A11 A12;   =  [R11 R12;
%       A21 A22]       0  R22], 
% where Q is orthogonal such that R11 (s by s) is upper triangular matrix.
% 
% This code allows for computing the post-array R and its derivative R'_{\theta} 
%
% References: (see Lemma 2 and Algorithm 2 from the paper)
%        Kulikova M.V., Tsyganova J.V. (2015) "Constructing numerically stable 
%        Kalman filter-based algorithms for gradient-based adaptive filtering", 
%        International Journal of Adaptive Control and Signal Processing, 
%        29(11):1411-1426. DOI http://dx.doi.org/10.1002/acs.2552 
%
% Input:
%     Pre_array      - given pre-array A;
%     s              - size of the block (R_11) that should be triangularized;
%     Diff_pre_array - given derivative of the pre-array A'_{\theta}; 
% Output:
%     Post_array      - get post-array R (block upper triangular);
%     Diff_post_array - get the derivative of the post-array R'_{\theta};
%

function [Post_array,Diff_post_array,Orthog] = Diff_QR(Pre_array,s,Diff_pre_array)

[sk,sl] = size(Pre_array); k = sk-s; l = sl-s;   % sizes of the corresponding blocks

[Orthog,Post_array] = qr(Pre_array(1:s+k,1:s));  % Q is any orthogonal rotation that
                                                 % (s+k) by (s) block of the pre-array                                                
if l~=0
 Post_array2 = Orthog'*Pre_array(1:end,end); % compute R12, R22 blocks
 Post_array  = [Post_array, Post_array2];    % the full post-array
end

Diff_post_array = [];
if nargin>2
   R11 = Post_array(1:s,1:s);                % blocks of the post-array
   R12 = Post_array(1:s,s+1:s+l);            % blocks of the post-array
   R22 = Post_array(s+1:s+k,s+1:s+l);        % blocks of the post-array 
      
   Q_applied_DiffA = Orthog'*Diff_pre_array; % apply Q to the pre-array derivatives
   XX = Q_applied_DiffA(1:s,1:s);            % notations
   YY = Q_applied_DiffA(s+1:s+k,1:s);
   NN = Q_applied_DiffA(1:s,s+1:s+l);
     
   matrix_product = XX/R11; 
   % compute the matrix product
   % Note, MatLab treats correctly the empty blocks
                         
   Diff_R11 = (tril(matrix_product,-1)'+triu(matrix_product))*R11;
   Diff_R12 = (tril(matrix_product,-1)'-tril(matrix_product,-1))*R12+R11'\YY'*R22+NN;
   Diff_R21 =  zeros(k,s);  
   Diff_R22 =  NaN(k,l); % The proposed method does not allow for 
                         % computing R22'_{\theta}. We set this block to NaN. 
                         % Fortunately, in the most filtering methods
                         % this block is of no interest
                         
   Diff_post_array = [Diff_R11, Diff_R12; Diff_R21, Diff_R22]; % Result
end
end
