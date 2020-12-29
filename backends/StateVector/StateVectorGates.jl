
module StateVectorGates

    export SV_Gate, SV_OneQubitGate, SV_TwoQubitGate
    abstract type SV_Gate end
    abstract type SV_OneQubitGate <: SV_Gate end
    abstract type SV_TwoQubitGate <: SV_Gate end

    struct X <: SV_OneQubitGate
        _target::Int64
    end

    struct Y <: SV_OneQubitGate
        _target::Int64
    end

    struct Z <: SV_OneQubitGate
        _target::Int64
    end

    struct H <: SV_OneQubitGate
        _target::Int64
    end

    struct T <: SV_OneQubitGate
        _target::Int64
    end

    struct S <: SV_OneQubitGate
        _target::Int64
    end

    struct RX <: SV_OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct RY <: SV_OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct RZ <: SV_OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct U3 <: SV_OneQubitGate
        _target::Int64
        _theta::Float64
        _phi::Float64
        _lambd::Float64
    end

    struct CZ <: SV_TwoQubitGate
        _control::Int64
        _target::Int64
    end

    struct CX <: SV_TwoQubitGate
        _control::Int64
        _target::Int64
    end

    struct CP <: SV_TwoQubitGate
        _control::Int64
        _target::Int64
    end

    struct CRX <: SV_TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
    end

    struct CRY <: SV_TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
    end

    struct CRZ <: SV_TwoQubitGate
        _control::Int64
        _target::Int64
        _theta::Float64
    end

    # export all SV_Gate and TwoQubitGate
    for n in names(@__MODULE__; all=true)
        if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
            if @eval typeof($n) <: DataType
                if @eval $n <: SV_Gate
                    @eval export $n
                end
            end
        end
    end
end