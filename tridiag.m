function y = tridiag( a, b, c, f )
%  This is an implementation of the Thomas algorithm.
%  %a*C(i-1) + b*C(i) + c*C(i+1) = f
%  Solve the  n x n  tridiagonal system for y:
%
%  [ b(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
%  [ a(2)  b(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
%  [       a(3)  b(3)  c(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        a(n-1) b(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
%  [                                 a(n)  b(n)  ] [  y(n)  ]   [  f(n)  ]
%
%  f must be a vector (row or column) of length n
%  a, b, c must be vectors of length n (note that a(1) and c(n) are not used)

n = length(f);
v = zeros(n,1);   
y = v;
w = b(1);
y(1) = f(1)/w;
for i=2:n
    v(i-1) = c(i-1)/w;
    w = b(i) - a(i)*v(i-1);
    y(i) = ( f(i) - a(i)*y(i-1) )/w;
end
for j=n-1:-1:1
   y(j) = y(j) - v(j)*y(j+1);
end



