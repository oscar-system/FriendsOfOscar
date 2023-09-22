using Oscar
using HomotopyContinuation
using Gapjm
using polyhedral_jll

# To translate OSCAR polynomials to something HomotopyContinuation can use
# courtesy of thofma
function to_dynamic_polynomial(f, vars = (@polyvar x[1:ngens(parent(f))])[1])
  @assert coefficient_ring(f) === QQ
  n = ngens(parent(f))
  return sum(c * prod(vars[i]^Int(e[i]) for i in 1:n) for (c, e) in zip(Oscar.coefficients(f), exponents(f)))
end
