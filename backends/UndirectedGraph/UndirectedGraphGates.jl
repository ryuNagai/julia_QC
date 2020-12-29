
module UndirectedGraphGates

    export OneQubitGateTensor, TwoQubitGateTensor, UG_Gate, UG_OneQubitGate, UG_TwoQubitGate, UndirectedGraphTensor
    abstract type UndirectedGraphTensor end
    abstract type UG_Gate end
    abstract type UG_OneQubitGate <: UG_Gate end
    abstract type UG_TwoQubitGate <: UG_Gate end

    mutable struct OneQubitGateTensor <: UndirectedGraphTensor
        item::Array{Complex, 2}
        target::Int
        worldline1::Int
        worldline2::Int
    end

    mutable struct TwoQubitGateTensor <: UndirectedGraphTensor
        item::Array{Complex, 2}
        control::Int
        target::Int
        worldline1::Int
        worldline2::Int
    end

    mutable struct X <: UG_OneQubitGate
        _target::Int64
    end

    mutable struct Y <: UG_OneQubitGate
        _target::Int64
    end

    mutable struct Z <: UG_OneQubitGate
        _target::Int64
    end

    mutable struct H <: UG_OneQubitGate
        _target::Int64
    end

    mutable struct T <: UG_OneQubitGate
        _target::Int64
    end

    struct S <: UG_OneQubitGate
        _target::Int64
    end

    struct RX <: UG_OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct RY <: UG_OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct RZ <: UG_OneQubitGate
        _target::Int64
        _theta::Float64
    end

    struct U3 <: UG_OneQubitGate
        _target::Int64
        _theta::Float64
        _phi::Float64
        _lambd::Float64
    end

    struct CZ <: UG_TwoQubitGate
        _control::Int64
        _target::Int64
    end

    struct CX <: UG_TwoQubitGate
        _control::Int64
        _target::Int64
    end

    struct CP <: UG_TwoQubitGate
        _control::Int64
        _target::Int64
    end

    # export all UG_Gate and TwoQubitGate
    for n in names(@__MODULE__; all=true)
        if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
            if @eval typeof($n) <: DataType
                if @eval $n <: UG_Gate
                    @eval export $n
                end
            end
        end
    end
end