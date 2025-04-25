module Honeycomb

using Oscar
# found at https://github.com/taboege/OscarHomotopyContinuation
using OscarHomotopyContinuation
using HomotopyContinuation
using Distributions

export honeycomb_meta_system, rand_honeycomb_system, record_data, record_data_numeric, record_data_solve

# lifting functions 
function lift_func(i, j)
    return i^2 + j^2 + i*j
end

function lift_func(l)
    i, j = l[1], l[2]
    return lift_func(i, j)
end

function honeycomb_meta_system(delta, l_func = "min")
    R, (x, y, t) = polynomial_ring(QQ, ["x", "y", "t"])

    P = delta * simplex(2)
    L = lattice_points(P)
    init_weights = (lift_func).(L)
    subdiv = subdivision_of_points(L, init_weights)

    weights = l_func == "min" ? min_weights(subdiv) : init_weights

    F = zeros(R, 3)

    for (l, w) in zip(L, weights)

        i, j = Int(l[1]), Int(l[2])

        term = t^w*x^i*y^j
        
        F[mod(i-j, 3) + 1] += term
    end

    return Oscar.ideal(R, F)
end

function rand_coefficients(k, h, d, targ_ring)

    x = rand(Uniform(k, h), d)
    return [targ_ring(convert(Rational, num)) for num in x]
end

function rand_honeycomb_system(delta, l_func,
                               t_bot, t_top, c_bot, c_top)

    I = honeycomb_meta_system(delta, l_func)
    R = base_ring(I)
    S, (x, y) = polynomial_ring(QQ, ["x", "y"])
    phi = hom(R, S, [x, y, rand_coefficients(t_bot, t_top, 1, S)...])

    c1 = rand_coefficients(c_bot, c_top, 3, S)
    c2 = rand_coefficients(c_bot, c_top, 3, S)

    M = gens(I)
    f1 = sum(c1 .* (phi).(M))
    f2 = sum(c2 .* (phi).(M))

    return ideal(S, [f1, f2])
end

# assumes that the deformation parameter of system is the last variable
# computes the elimnant in this parameter and returns rational approx.
# of its real roots
function elim_meta_system(meta_sys_ideal)

    R = base_ring(meta_sys_ideal)
    n = ngens(R)
    gb = groebner_basis_f4(meta_sys_ideal, eliminate = n - 1,
                           complete_reduction = true,
                           info_level = 2)

    elimn = first(gb)
    elimn_sys = Oscar.ideal(R, [elimn, gens(R)[1:n-1]...])
    rsols = Oscar.real_solutions(elimn_sys)

    return (elimn, (x -> last(x)).(rsols[1]))
end

function record_data(delta_max, l_func)

    for delta in 3:delta_max
        iseven(delta) && continue
        f = open("./msolve_elim_$(delta)_$(l_func).txt", "w+")
        println("delta = $(delta)")

        I = honeycomb_meta_system(delta, l_func)
        try
            tim = @elapsed elimn, rroots = elim_meta_system(I)
            println(f, "delta = $(delta)")
            println(f, "timing $(tim)")
            println(f, elimn)
            println(f, "$(length(rroots)) real roots")
            println(f, "real roots from elimnant:")
            println(f, rroots)
            
            println(f, "------")
            close(f)
        catch e
            println(f, e)
            close(f)
            continue
        end

    end
end

function record_data_solve(delta_max, l_func)

    for delta in 3:delta_max
        iseven(delta) && continue
        f = open("./msolve_solve_$(delta)_$(l_func).txt", "w+")
        println("delta = $(delta)")
        I = honeycomb_meta_system(delta, l_func)

        tim2 = 0.0
        try
            tim2 = @elapsed rroots2 = (x -> last(x)).(Oscar.real_solutions(I, info_level = 2)[1])

            println(f, "delta = $(delta)")
            println(f, "timing from solving: $(tim2)")
            println(f, "$(length(rroots2)) real roots")
            println(f, "real roots from solving:")
            println(f, rroots2)
            
            println(f, "------")
            close(f)
        catch e
            println(f, e)
            close(f)
            continue
        end

    end 
end

function record_data_numeric(delta_max, l_func)

    for delta in 3:delta_max
        iseven(delta) && continue
        f = open("./hc_solve_$(delta)_$(l_func).txt", "w+")
        println("delta = $(delta)")

        I = honeycomb_meta_system(delta, l_func)
        sys = OscarHomotopyContinuation.System(I)

        tim1 = tim2 = 0.0
        try
            tim1 = @elapsed sls = HomotopyContinuation.solve(sys, show_progress = true,
                                                             start_system = :polyhedral,
                                                             threading = false)
            tim2 = @elapsed certs = certify(sys, sls, threading = false)
            tim = tim1 + tim2

            n_total = nsolutions(sls)
            
            n_real = nreal(sls)
            n_real_singular = nsingular(sls, only_real = true)
            
            n_real_cert = nreal_certified(certs)
            real_certs = filter(HomotopyContinuation.is_real, certificates(certs))
            
            real_certs_ints = (certified_solution_interval).(real_certs)
            real_certs_t_intervall = (x -> x[1,1]).(real_certs_ints)
            
            println(f, "delta = $(delta)")
            println(f, "timing $(tim)")
            println(f, "$(n_real_cert)/$(n_real) certified real solutions")
            println(f, "$(n_real_singular) singular real solutions")
            println(f, real_certs_t_intervall)
            println(f, "------")
            close(f)
        catch e
            println(f, e)
            close(f)
            continue
        end

    end
end

end # module Honeycomb
