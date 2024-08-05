using OscarHomotopyContinuation
using Oscar

R, (x,y) = QQ["x", "y"]
I = ideal(R, [ 2*x^2 + y^2 - 1, x^2 + 2*y^2 - 1 ])
@show nsolve(I)
@show ndim(I)
@show ndegree(I)
