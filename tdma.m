
function x = tdma(A,d)
% tdma: Tridiaogonal Matrix solver
% x = tdma(A,d):
%   Solves a tridiagonal matrix using the Thomas Algorithm
%
% input: 
%   A = Coefficient matrix (A square tridiagonal matrix)
%   d = Constant vector
% output:
%   x = Solution vector
%
% Author: Divyaprakash, PhD Candidate
%         Applied Mechanics Department, IIT Delhi
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 16 October 2020
    
% Number of unknowns
N = numel(d);
% Initialize the solution vector
x = zeros(N,1);
% Extract the diagonals and make their length equal to N by padding with 0's
a = [0; diag(A,-1)];
b = diag(A);
c = [diag(A,1); 0];
% Thomas Algorithm
for i = 1:N
    if i == 1
        c(i) = c(i)/b(i);
        d(i) = d(i)/b(i);  
    else
        c(i) = c(i)/(b(i) - a(i)*c(i-1));
        d(i) = (d(i) - a(i)*d(i-1))/(b(i) - a(i)*c(i-1));
    end
end
for i = N:-1:1
    if i == N
        x(i) = d(i);
    else
        x(i) = d(i) - c(i)*x(i+1); 
    end
end
end
