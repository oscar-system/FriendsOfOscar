# utility function:
# character values are represented by `QQAbFieldElem`s;
# return the corresponding `BigInt` for integral values,
# return the corresponding `QQBarFieldElem` for irrational values
function _map_irrat_to_field(F::QQBarField, x::QQAbFieldElem)
  N = x.c
  data = x.data
  cfs = coefficients(data)
  N == 1 && return BigInt(ZZ(cfs[1]))
  r = root_of_unity(F, N)
  pow = one(F)
  res = cfs[1] * pow
  for i in 2:length(cfs)
    pow = pow * r
    res = res + cfs[i] * pow
  end
  return res
end


# utility function:
# For a vector `A` of vectors (class functions),
# compute a vector for which duplicate columns are removed,
# and return `M, fus, proj`
# where `M` is the vector of collapsed vectors,
# and `fus` and `proj` are vectors of integers
# that describe the relation between `A` and `M`.
function _collapse_columns(A::Vector)
  m = length(A)
  n = length(A[1])
  fus = collect(1:n)
  proj = Int[1]
  take = ones(Bool, n)
  for j in 2:n
    # find the first column equal to the j-th column
    for k in 1:(j-1)
      if all(v -> v[k] == v[j], A)
        fus[j] = fus[k]
        take[j] = false
        continue
      end
    end
    if take[j]
      push!(proj, j)
      fus[j] = length(proj)
    end
  end

  return map(v -> v[findall(take)], A), fus, proj
end


"""
    s_character_simplex(tbl::Oscar.GAPGroupCharacterTable;
               rational::Bool = true,
               ppow_nonzero::Bool = false,
               irrats::Symbol = :nf)

Return `P, galoissums, ppow_pos` where
- `P` is the polyhedron that is defined by the inequalities given by
  the rational irreducible characters of `tbl` (if `rational` is `true`)
  or the real irreducible characters of `tbl` (if `rational` is `false`),

- `galoissums` is the vector of real or rational irreducibles of `tbl`, and

- `ppow_pos` is a vector of column positions in `tbl` for which a strictly
  positive value is prescribed in the S-characters of `tbl` encoded by `P`;
  these are the rational classes of elements of prime power order if
  `ppow_nonzero` is `true`, otherwise `ppow_pos` is empty.

The optional argument `irrats` defines how the irrational entries of the
defining matrix of `P` are represented:
The value `:nf` means that these entries belong to a common embedded
number field (constructed from the common field of character values),
other values mean that these entries lie in `algebraic_closure(QQ)`.

The polytope `P` is known to be a simplex.
"""
function s_character_simplex(tbl::Oscar.GAPGroupCharacterTable;
                    rational::Bool = true,
                    ppow_nonzero::Bool = false,
                    irrats::Symbol = :nf)
  @req order(tbl) > 1 "the underlying group must not be trivial"

  # Compute the orbit sums on the nontrivial irreducibles
  # under complex conjugation if `rational == false`
  # and under all Galois automorphisms if `rational == true`;
  # omit the trivial character.
  triv = trivial_character(tbl)
  pos = findfirst(isequal(triv), tbl)
  if rational
    c = collect(tbl)
    nontriv = vcat(c[1:pos-1], c[pos+1:end])
    galoissums = unique(map(galois_orbit_sum, nontriv))
  else
    galoissums = Oscar.GAPGroupClassFunction[]
    for i in 1:length(tbl)
      i == pos && continue
      chi = tbl[i]
      psi = conj(chi)
      if psi == chi
        push!(galoissums, chi)
      elseif findfirst(isequal(psi), tbl) > i
        push!(galoissums, chi + psi)
      end
    end
  end

  # Collapse equal columns (do not duplicate the defining inequalities).
  coll, fus, proj = _collapse_columns(map(values, galoissums))

  # Negate and transpose the matrix,
  # map irrational entries to the field in question.
  m = length(proj)
  n = length(coll)
  A = Matrix(undef, m, n)
  ratcol = ones(Bool, m)
  if irrats == :nf
    CF, emb = character_field(galoissums)
    e = real_embeddings(CF)[1]
    F = Hecke.embedded_field(CF, e)[1]
    for i in 1:m, j in 1:n
      x = coll[j][i]
      if x.c == 1
        A[i,j] = - BigInt(ZZ(coefficients(x.data)[1]))
      else
        ratcol[i] = false
        A[i,j] = - F(preimage(emb, x))
      end
    end
  else
    F = algebraic_closure(QQ)
    for i in 1:m, j in 1:n
      val = coll[j][i]
      if val.c != 1
        ratcol[i] = false
      end
      A[i,j] = - _map_irrat_to_field(F, val)
    end
  end

  # Initialize the right hand side `b` of the system of inequalities.
  b = ones(Int, m)
  if ppow_nonzero
    # Replace the `1` in `b` by `0` in all those positions of columns
    # that are rational and belong to elements of prime power order.
    ppow_pos = findall(x -> x == 1 || is_prime_power_with_data(x)[1],
                       orders_class_representatives(tbl))
    for i in unique(fus[ppow_pos])
      if rational || ratcol[i]
        b[i] = 0
      end
    end
  else
    # Do not prescribe a stronger condition for elements of prime power order.
    ppow_pos = Int[]
  end

  # Compute the polyhedron that is defined by the inequalities
  # given by the (real or rational) irreducible characters.
  # Its lattice points are the coefficient vectors of
  # all integer linear combinations `psi` of the real or rational irreducibles
  # such that `psi + 1` has nonnegative values
  # (and if applicable then `psi + 1` is nonzero on certain classes,
  # see above).
  P = polyhedron(F, A, b)

  return P, galoissums, ppow_pos
end


"""
    s_characters(tbl::Oscar.GAPGroupCharacterTable;
                 rational::Bool = true,
                 ppow_nonzero::Bool = false,
                 irrats::Symbol = :nf)

Return the vector of nontrivial S-characters of `tbl`
that are described by the arguments.
See [`s_character_simplex`](@ref) for the meaning of the keyword arguments.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> res = s_characters(tbl, rational = false, ppow_nonzero = false);

julia> length(res)  # all nontrivial S-characters of A5
24

julia> res = s_characters(tbl, rational = true, ppow_nonzero = false);

julia> length(res)  # all rational nontrivial S-characters of A5
16

julia> res = s_characters(tbl, rational = false, ppow_nonzero = true);

julia> length(res)  # no nontrivial S-characters positive on p-power elements
0

julia> res = s_characters(tbl, rational = true, ppow_nonzero = true);

julia> length(res)  # of course no rational such S-characters
0
```
"""
function s_characters(tbl::Oscar.GAPGroupCharacterTable;
                       rational::Bool = true,
                       ppow_nonzero::Bool = false,
                       irrats::Symbol = :nf)
  # The trivial group has no nontrivial S-characters.
  order(tbl) == 1 && return typeof(tbl[1])[]

  P, galoissums, ppow_pos = s_character_simplex(tbl,
                              rational = rational,
                              ppow_nonzero = ppow_nonzero,
                              irrats = irrats)
  ll = lattice_points(P)

  # Compute the corresponding S-characters;
  # if `ppow_nonzero` is true then collect those that are nonzero
  # on all classes of elements of prime power order.
  res = typeof(galoissums[1])[]
  triv = trivial_character(parent(galoissums[1]))
  n = length(galoissums)
  for v in ll
    is_zero(v) && continue
    cand = triv + sum(v[i] * galoissums[i] for i in 1:n)
    if !ppow_nonzero || all(x -> !is_zero(cand[x]), ppow_pos)
      push!(res, cand)
    end
  end

  return res
end

function s_characters_via_factor_group(tbl::Oscar.GAPGroupCharacterTable,
             nclasses::Vector{Int};
             rational::Bool = true,
             ppow_nonzero::Bool = false,
             irrats::Symbol = :nf)

  # Compute the S-characters for the factor group in question.
  facttbl, factfus = quo(tbl, nclasses)
  fact_s_chars = s_characters(facttbl,
                       rational = rational,
                       ppow_nonzero = ppow_nonzero,
                       irrats = irrats)

  # Compute the orbit sums on the N-faithful irreducibles
  # under complex conjugation if `rational == false`
  # and under all Galois automorphisms if `rational == true`.
  if rational
    N_faith = filter(x -> !is_subset(nclasses, class_positions_of_kernel(x)),
                     collect(tbl))
    galoissums = unique(map(galois_orbit_sum, N_faith))
  else
    galoissums = Oscar.GAPGroupClassFunction[]
    for i in 1:length(tbl)
      chi = tbl[i]
      is_subset(nclasses, class_positions_of_kernel(chi)) && continue
      psi = conj(chi)
      if psi == chi
        push!(galoissums, chi)
      else
        pos = findfirst(isequal(psi), tbl)
        if pos > i
          push!(galoissums, chi + psi)
        end
      end
    end
  end

  # Collapse equal columns (do not duplicate the defining inequalities).
  coll, fus, proj = _collapse_columns(map(values, galoissums))

  # Negate and transpose the matrix,
  # map irrational entries to the field in question.
  m = length(proj)
  n = length(coll)
  A = Matrix(undef, m, n)
  ratcol = ones(Bool, m)
  if irrats == :nf
    # We may need a proper field extension.
    # Make sure that the field contains also
    # the character values for the factor group.
    # (Since we consider only *one* solution for the factor group
    # at a time, it is only necessary to consider the field of values
    # of this vector.
    # All solutions we found up to now are rational,
    # thus we need not think about such an improved strategy yet.)
    CF, emb = character_field(vcat(galoissums, fact_s_chars))
    e = real_embeddings(CF)[1]
    F = Hecke.embedded_field(CF, e)[1]
    for i in 1:m, j in 1:n
      x = coll[j][i]
      if x.c == 1
        A[i,j] = - BigInt(ZZ(coefficients(x.data)[1]))
      else
        ratcol[i] = false
        A[i,j] = - F(preimage(emb, x))
      end
    end
  else
    F = algebraic_closure(QQ)
    for i in 1:m, j in 1:n
      val = coll[j][i]
      if val.c != 1
        ratcol[i] = false
      end
      A[i,j] = - _map_irrat_to_field(F, val)
    end
  end

  if ppow_nonzero
    # Decrease the right hand side by `1` in all those positions of columns
    # that are rational and belong to elements of prime power order.
    ppow_pos = findall(x -> x == 1 || is_prime_power_with_data(x)[1],
                       orders_class_representatives(tbl))
  else
    # Do not prescribe a stronger condition for elements of prime power order.
    ppow_pos = Int[]
  end

  # For each solution v, compute the solutions in the next smaller
  # factor group that lie above v.
  res = [restrict(chi, tbl) for chi in fact_s_chars]

  # Solutions for `tbl` may project to the trivial character of the factor.
  pushfirst!(fact_s_chars, trivial_character(facttbl))

  for v in fact_s_chars
    # Initialize the right hand side `b` of the system of inequalities.
    # (Fill up v by characters not from the factor group)
    # First inflate v, then collapse it (always take the strongest condition!).
#T really do this!
    infl_v = values(v)[factfus]
    coll_v = infl_v[proj]
    if irrats == :nf
      b = map(x -> F(preimage(emb,x)), coll_v)
    else
      b = map(x -> _map_irrat_to_field(F, x), coll_v)
    end
    if ppow_nonzero
      for i in unique(fus[ppow_pos])
        if rational || ratcol[i]
          b[i] = b[i] - 1
        end
      end
    end

    P = polyhedron(F, A, b)
    ll = lattice_points(P)

    # Compute the corresponding S-characters;
    # if `ppow_nonzero` is true then collect those that are nonzero
    # on all classes of elements of prime power order.
    factcand = Oscar.class_function(tbl, infl_v)
    n = length(galoissums)
    for w in ll
      is_zero(w) && continue
      cand = factcand + sum(w[i] * galoissums[i] for i in 1:n)
      if !ppow_nonzero || all(x -> !is_zero(cand[x]), ppow_pos)
        push!(res, cand)
      end
    end
  end

  return res
end;


# utility function:
# save the result of `s_characters` to a file
function save_s_characters(filename::String, s_chars::Vector)
  l = [coordinates(ZZRingElem, v) for v in s_chars]
  tbl = parent(s_chars[1])
  info = Dict(:name => identifier(tbl), :coordinates => l)
  save(filename, info)
end


# utility function:
# create the LaTeX table in the paper from the stored data
names = [
          [ # alternating groups and their extensions
            [ "A8", "\\fA_{8}" ],
            [ "2.A8", "2.\\fA_{8}" ],
            [ "A9", "\\fA_{9}" ],
            [ "2.A9", "2.\\fA_{9}" ],
            [ "A10", "\\fA_{10}" ],
            [ "2.A10", "2.\\fA_{10}" ],
            [ "A11", "\\fA_{11}" ],
            [ "2.A11", "2.\\fA_{11}" ],
          ],
          [ # sporadic simple groups and their extensions
            [ "M12", "M_{12}" ],
            [ "2.M12", "2.M_{12}" ],
            [ "M12.2", "M_{12}.2" ],
            [ "J1", "J_{1}" ],
            [ "J2", "J_{2}" ],
            [ "2.J2", "2.J_{2}" ],
            [ "J2.2", "J_{2}.2" ],
            [ "J3", "J_{3}" ],
            [ "McL", "M^{c}L" ],
            [ "HS", "HS" ],
            [ "2.HS", "2.HS" ],
            [ "M24", "M_{24}" ],
          ],
          [ # simple groups of Lie type and their extensions
            [ "L4(3)", "\\PSL_{4}(3)" ],
            [ "L5(2)", "\\PSL_{5}(2)" ],
            [ "L5(2).2", "\\PSL_{5}(2).2" ],
            [ "S4(4)", "\\PSp_{4}(4)" ],
            [ "S4(4).2", "\\PSp_{4}(4).2" ],
            [ "S4(5)", "\\PSp_{4}(5)" ],
            [ "S6(2)", "\\PSp_{6}(2)" ],
            [ "U4(2)", "\\PSU_{4}(2)" ],
            [ "2.U4(2)", "2.\\PSU_{4}(2)" ],
            [ "2F4(2)'", "{}^{2}F_{4}(2)'" ],
            [ "2F4(2)'.2", "{}^{2}F_{4}(2)'.2" ],
            [ "R(27)", "^2G_2(27)" ],
            [ "G2(3)", "G_{2}(3)" ],
          ],
          [ # some other perfect groups
            [ "2^4:a8", "2^{4}\\!:\\!\\fA_{8}" ],
            [ "2^4.a8", "2^{4}.\\fA_{8}" ],
            [ "2^5.psl(5,2)", "2^{5}.\\PSL_{5}(2)" ],
            [ "2^6:A8", "2^{6}\\!:\\!\\fA_{8}" ],
          ],
        ];

function s_characters_table(names)
  # table header
  println("\\begin{table}")
  println("\\caption{Examples}   \\label{tab:ex}")
  println("\\(\n\\begin{array}{lrrrrr}\n  \\toprule")
  println("  G & \\# \\text{classes} & \\# \\text{real} & \\# \\text{rat.} & ",
          "\\# \\text{S-char.}\n",
          "  & \\# \\text{virt.\\ S-char.} \\\\\n",
          "  \\midrule")

  # run over the name portions
  for l in names
    for pair in l
      name = pair[1]
      latexname = pair[2]

      r = load(joinpath(@__DIR__, "../data", name))
      list = r[:coordinates]
      tbl = character_table(name)

      # compute the *faithful* example characters,
      # i.e. those which do not live in a proper factor group.
      # Note that we cannot use `class_positions_of_kernel` for
      # virtual characters.
      nsg = setdiff(class_positions_of_normal_subgroups(tbl), [[1]])
      nccl = length(tbl)
      faithirrpos = [filter(i -> !is_subset(n, class_positions_of_kernel(tbl[i])), 1:nccl) for n in nsg]
      faith = filter(x -> all(l -> any(i -> x[i] != 0, l), faithirrpos), list)

      # compute the *faithful* ordinary ones:
      ord = filter(x -> minimum(x) >= 0, list)
      ordfaith = intersect(ord, faith)

      # compute the values for the first three columns
      nreal = length(unique([x + conj(x) for x in tbl]))
      nrat = length(unique(map(galois_orbit_sum, collect(tbl))))

      # print the table row for this group,
      # separate portions
      if length(faith) == 0
        print("% ")
      else
        print("  ")
      end
      println(latexname, " & ", nccl, " & ", nreal, " & ", nrat, " & ",
              length(faith), " & ", length(faith) - length(ordfaith), " \\\\",
              (pair == l[end] && l != names[end]) ? "%\n  [1.6ex]\n  %" : "" )
    end
  end

  # table footer
  println("  \\bottomrule\n\\end{array}\n\\)\n\\end{table}")
end
