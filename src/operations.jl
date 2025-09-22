function removenode!(n::Node, t::Node, h::Node, v::Vehicle, s::Solution)
    # fetch details
    x = v.x * v.n
    y = v.y * v.n
    aₜₙ = s.A[t.i, n.i]
    aₙₕ = s.A[n.i, h.i]
    aₜₕ = s.A[t.i, h.i]
    δ = (aₜₙ.c + aₙₕ.c) - aₜₕ.c
    # penalty
    s.p -= (v.l > v.q) ? (v.l - v.q) : 0
    # update node
    n.t = 0
    n.h = 0
    n.v = 0
    t.h = isdepot(t) ? t.h : h.i
    h.t = isdepot(h) ? h.t : t.i
    # update vehicle
    v.s = isdepot(t) ? h.i : v.s
    v.e = isdepot(h) ? t.i : v.e
    v.n -= 1
    v.l -= n.q
    v.x = iszero(v.n) ? 0. : (x - n.x) / v.n
    v.y = iszero(v.n) ? 0. : (y - n.y) / v.n
    v.c -= δ
    # update solution
    s.c -= δ
    #penalty
    s.p += (v.l > v.q) ? (v.l - v.q) : 0
    # return solution
    return s
end

function insertnode!(n::Node, t::Node, h::Node, v::Vehicle, s::Solution)
    x = v.x * v.n
    y = v.y * v.n
    aₜₙ = s.A[t.i, n.i]
    aₙₕ = s.A[n.i, h.i]
    aₜₕ = s.A[t.i, h.i]
    δ = (aₜₙ.c + aₙₕ.c) - aₜₕ.c
    # penalty
    s.p -= (v.l > v.q) ? (v.l - v.q) : 0
    # update node
    n.t = t.i
    n.h = h.i
    n.v = v.i
    t.h = isdepot(t) ? t.h : n.i
    h.t = isdepot(h) ? h.t : n.i
    # update vehicle
    v.s = isdepot(t) ? n.i : v.s
    v.e = isdepot(h) ? n.i : v.e
    v.n += 1
    v.l += n.q
    v.x = (x + n.x) / v.n
    v.y = (y + n.y) / v.n
    v.c += δ
    # update solution
    s.c += δ
    # penalty
    s.p += (v.l > v.q) ? (v.l - v.q) : 0
    # return solution
    
    return s
end