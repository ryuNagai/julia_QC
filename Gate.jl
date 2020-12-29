module GateSet
    export Gate

    abstract type Gate end
    abstract type OneQubitGate <: Gate end
    abstract type TwoQubitGate <: Gate end

    struct X <: OneQubitGate
        _target::Int64
        _name::String
    end

    function X(_target::Int64)
        return X(_target, "X")
    end

    struct Y <: OneQubitGate
        _target::Int64
        _name::String
    end

    function Y(_target::Int64)
        return Y(_target, "Y")
    end

    struct Z <: OneQubitGate
        _target::Int64
        _name::String
    end

    function Z(_target::Int64)
        return Z(_target, "Z")
    end

    struct H <: OneQubitGate
        _target::Int64
        _name::String
    end

    function H(_target::Int64)
        return H(_target, "H")
    end

    struct T <: OneQubitGate
        _target::Int64
        _name::String
    end

    function T(_target::Int64)
        return T(_target, "T")
    end

    struct S <: OneQubitGate
        _target::Int64
        _name::String
    end

    function S(_target::Int64)
        return S(_target, "S")
    end

    struct RX <: OneQubitGate
        _target::Int64
        _theta::Float64
        _name::String
    end

    function RX(_target::Int64, _theta::Float64)
        return RX(_target, _theta, "RX")
    end

    struct RY <: OneQubitGate
        _target::Int64
        _theta::Float64
        _name::String
    end

    function RY(_target::Int64, _theta::Float64)
        return RY(_target, _theta, "RY")
    end

    struct RZ <: OneQubitGate
        _target::Int64
        _theta::Float64
        _name::String
    end

    function RZ(_target::Int64, _theta::Float64)
        return RZ(_target, _theta, "RZ")
    end

    struct U3 <: OneQubitGate
        _target::Int64
        _theta::Float64
        _phi::Float64
        _lambd::Float64
        _name::String
    end

    function U3(_target::Int64, _theta::Float64, _phi::Float64, _lambd::Float64)
        return U3(_target, _theta, _phi, _lambd, "U3")
    end

    struct CX <: TwoQubitGate
        _control::Int64
        _target::Int64
        _name::String
    end

    function CX(_control::Int64, _target::Int64)
        return CX(_control, _target, "CX")
    end

    struct CZ <: TwoQubitGate
        _control::Int64
        _target::Int64
        _name::String
    end

    function CZ(_control::Int64, _target::Int64)
        return CZ(_control, _target, "CZ")
    end

    struct CP <: TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
        _name::String
    end

    function CP(_control::Int64, _target::Int64, _theta::Float64)
        return CP(_control, _target, _theta, "CP")
    end

    struct CRX <: TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
        _name::String
    end

    function CRX(_control::Int64, _target::Int64, _theta::Float64)
        return CRX(_control, _target, _theta, "CRX")
    end

    struct CRY <: TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
        _name::String
    end

    function CRY(_control::Int64, _target::Int64, _theta::Float64)
        return CRY(_control, _target, _theta, "CRY")
    end

    struct CRZ <: TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
        _name::String
    end

    function CRZ(_control::Int64, _target::Int64, _theta::Float64)
        return CRZ(_control, _target, _theta, "CRZ")
    end

    # export all OneQubit Gate and TwoQubitGate
    for n in names(@__MODULE__; all=true)
        if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
            if @eval typeof($n) <: DataType
                if @eval $n <: Gate
                    @eval export $n
                end
            end
        end
    end

end