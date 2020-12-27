module StateVectorFunctions
    include("../../BaseModules.jl")
    include("../../Gate.jl")

    using Reexport
    @reexport using .GateSet
    @reexport using .BaseModule
    export apply!

    """
    X gate
    """
    function apply!(gate::X, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i )
            t = qubits[i0]
            qubits[i0] = qubits[i0 + (1 << (target))]
            qubits[i0 + (1 << (target))] = t
        end
        return qubits
    end

    """
    Y gate
    """
    function apply!(gate::Y, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i)
            t = qubits[i0] * 1im
            qubits[i0] = qubits[i0 + (1 << target)] * -1im
            qubits[i0 + (1 << target)] = t
        end
        return qubits
    end

    """
    Z gate
    """
    function apply!(gate::Z, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            qubits[_shifted(lower_mask, i) + (1 << target)] *= -1
        end
        return qubits
    end

    """
    H gate
    """
    function apply!(gate::H, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        sqrt2_inv = 1 / sqrt(2)#0.7071067811865475
        lower_mask = (1 << target) - 1
        for i in prange(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i)
            t = qubits[i0]
            u = qubits[i0 + (1 << target)]
            qubits[i0] = (t + u) * sqrt2_inv
            qubits[i0 + (1 << target)] = (t - u) * sqrt2_inv
        end
        return qubits
    end

    """
    RX gate
    """
    function apply!(gate::RX, qubits::Array{Complex, 1}, n_qubits::Int)
        target = gate._target
        target = n_qubits - target
        ang = gate._theta * 0.5
        angpi = ang / pi
        _cos = cospi(angpi)
        nisin = sinpi(angpi) * -1im
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i)
            t = qubits[i0]
            u = qubits[i0 + (1 << target)]
            qubits[i0] = _cos * t + nisin * u
            qubits[i0 + (1 << target)] = nisin * t + _cos * u
        end
        return qubits
    end

    """
    RZ gate
    """
    function apply!(gate::RZ, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        ang = gate._theta * 0.5
        eit = exp(1im * ang)
        eitstar = conj(eit)
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i)
            t = qubits[i0]
            u = qubits[i0 + (1 << target)]
            qubits[i0] *= eitstar
            qubits[i0 + (1 << target)] *= eit
        end
        return qubits
    end

    function apply!(gate::RY, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        ang = gate._theta * 0.5
        angpi = ang / pi
        _cos = cospi(angpi)
        _sin = sinpi(angpi)
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i)
            t = qubits[i0]
            u = qubits[i0 + (1 << target)]
            qubits[i0] = _cos * t - _sin * u
            qubits[i0 + (1 << target)] = _sin * t + _cos * u
        end
        return qubits
    end

    function apply!(gate::U3, qubits::Array{Complex, 1}, n_qubits::Int)
        target = n_qubits - gate._target
        theta = gate._theta * 0.5
        thetapi = theta / pi
        phipi = gate._phi / pi
        lambdpi = gate._lambd / pi
        _cos = cos(thetapi)
        _sin = sin(thetapi)
        expadd = exp((phipi + lambdpi) * 0.5im)
        expsub = exp((phipi - lambdpi) * 0.5im)
        a = conj(expadd) * _cos
        b = conj(-expsub) * _sin
        c = expsub * _sin
        d = expadd * _cos
        lower_mask = (1 << target) - 1
        for i in 1:(1 << (n_qubits - 1))
            i0 = _shifted(lower_mask, i)
            t = qubits[i0]
            u = qubits[i0 + (1 << target)]
            qubits[i0] = a * t + b * u
            qubits[i0 + (1 << target)] = c * t + d * u
        end
        return qubits
    end
    
    function apply!(gate::CX, qubits::Array{Complex, 1}, n_qubits::Int)
        control = n_qubits - gate._control
        target = n_qubits - gate._target
        c_mask = 0
        c_mask |= 1 << control
        t_mask = 1 << target
        n_loop = 1 << (n_qubits - 2)
        masks = _create_masks([control, target])
        for i in 1:n_loop
            i10 = _mult_shifted(masks, i-1) | c_mask
            i11 = i10 | t_mask
            i10 += 1
            i11 += 1
            t = qubits[i10]
            qubits[i10] = qubits[i11]
            qubits[i11] = t
        end
        return qubits
    end

    function apply!(gate::CZ, qubits::Array{Complex, 1}, n_qubits::Int)
        control = n_qubits - gate._control
        target = n_qubits - gate._target
        all1 = 0
        for b in [control, target]
            all1 |= 1 << b
        end
        n_loop = 1 << (n_qubits - length([control, target]))
        masks = _create_masks([control, target])
        for i in 1:n_loop
            i11 = _mult_shifted(masks, i-1) | all1
            qubits[i11 + 1] *= -1
        end
        return qubits
    end

    function apply!(gate::CRX, qubits::Array{Complex, 1}, n_qubits::Int)
        control = n_qubits - gate._control
        target = n_qubits - gate._target
        ang = gate._theta * 0.5
        angpi = ang / pi
        _cos = cospi(angpi)
        nisin = sinpi(angpi) * -1.0im
        c_mask = 0
        c_mask |= 1 << control
        t_mask = 1 << target
        n_loop = 1 << (n_qubits - length([control, target]))
        masks = _create_masks([control, target])
        for i in 1:n_loop
            i10 = _mult_shifted(masks, i-1) | c_mask
            i11 = i10 | t_mask
            i10 += 1
            i11 += 1
            t = qubits[i10]
            u = qubits[i11]
            qubits[i10] = _cos * t + nisin * u
            qubits[i11] = nisin * t + _cos * u
        end
        return qubits
    end

    function apply!(gate::CRY, qubits::Array{Complex, 1}, n_qubits::Int)
        control = n_qubits - gate._control
        target = n_qubits - gate._target
        ang = gate._theta * 0.5
        angpi = ang / pi
        _cos = cospi(angpi)
        _sin = sinpi(angpi)
        c_mask = 0
        c_mask |= 1 << control
        t_mask = 1 << target
        n_loop = 1 << (n_qubits - length([control, target]))
        masks = _create_masks([control, target])
        for i in 1:n_loop
            i10 = _mult_shifted(masks, i-1) | c_mask
            i11 = i10 | t_mask
            i10 += 1
            i11 += 1
            t = qubits[i10]
            u = qubits[i11]
            qubits[i10] = _cos * t - _sin * u
            qubits[i11] = _sin * t + _cos * u
        end
        return qubits
    end

    function apply!(gate::CRZ, qubits::Array{Complex, 1}, n_qubits::Int)
        control = n_qubits - gate._control
        target = n_qubits - gate._target
        ang = gate._theta * 0.5
        eit = exp(1.0im * ang)
        eitstar = conj(eit)
        c_mask = 0
        c_mask |= 1 << control
        t_mask = 1 << target
        n_loop = 1 << (n_qubits - length([control, target]))
        masks = _create_masks([control, target])
        for i in 1:n_loop
            i10 = _mult_shifted(masks, i-1) | c_mask
            i11 = i10 | t_mask
            i10 += 1
            i11 += 1
            qubits[i10] *= eitstar
            qubits[i11] *= eit
        end
        return qubits
    end

    function apply!(gate::CP, qubits::Array{Complex, 1}, n_qubits::Int)
        control = n_qubits - gate._control
        target = n_qubits - gate._target
        eit = exp(1.0im * gate._theta)
        c_mask = 0
        c_mask |= 1 << control
        t_mask = 1 << target
        n_loop = 1 << (n_qubits - length([control, target]))
        masks = _create_masks([control, target])
        for i in 1:n_loop
            i11 = _mult_shifted(masks, i-1) | c_mask | t_mask + 1
            qubits[i11] *= eit
        end
        return qubits
    end

    function _create_masks(indices)
        sort!(indices)
        masks = zeros(Int, length(indices) + 1)
        for (i, x) in enumerate(indices)
            masks[i] = (1 << (x - (i - 1))) - 1
        end
        masks[length(masks)] = ~0
        for i in reverse(1:length(indices))
            masks[i + 1] &= ~masks[i]
        end
        return masks
    end
    
    function _mult_shifted(masks, idx)
        shifted = 0
        for (i, x) in enumerate(masks)
            shifted |= (idx & x) << (i-1)
        end
        return shifted
    end
end