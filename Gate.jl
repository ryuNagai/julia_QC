module GateSet
    export Gate, OneQubitGate, TwoQubitGate

    abstract type Gate end
    abstract type OneQubitGate <: Gate end
    abstract type TwoQubitGate <: Gate end


    struct X <: OneQubitGate
        _target::Int64
    end

    struct Y <: OneQubitGate
        _target::Int64
    end

    struct Z <: OneQubitGate
        _target::Int64
    end

    struct H <: OneQubitGate
        _target::Int64
    end

    struct T <: OneQubitGate
        _target::Int64
    end

    struct S <: OneQubitGate
        _target::Int64
    end

    struct RX <: OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct RY <: OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct RZ <: OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct CX <: TwoQubitGate
        _control::Int64
        _target::Int64
    end

    struct CZ <: TwoQubitGate
        _control::Int64
        _target::Int64
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