module Juliaqat

    #include("./Gate.jl")
    #include("./BaseModules.jl")
    include("./backends/backend_functions.jl")
    using Reexport
    @reexport using .BackendFunctions

    export QuantumCircuit, Input, apply!, get_device, execute

    devices = Dict([("UndirectedGraph", UndirectedGraphModel),
                    ("StateVector", StateVectorModel)
                        ])

    mutable struct QuantumCircuit
        _n_qreg::Int64
        _n_creg::Int64
        _gates::Array{Gate,1}
    end

    mutable struct Input
        _n_qreg::Int64
        _init_state::Array{Complex, 1}
        _gates::Array{Gate,1}
    end

    function get_device(name::String, x...)
        return devices[name](x)
    end

    function QuantumCircuit(n_qubit::Int64, n_creg::Int64)
        _n_qubit = n_qubit
        _n_creg = n_creg
        return QuantumCircuit(_n_qubit, _n_creg, [])
    end

    function Input(init_state::Array{Complex, 1})
        _init_state = init_state
        _n_qreg = length(init_state)
        return Input(n_qreg, init_state, [])
    end

    function Input(n_qreg::Int)
        _init_state = zeros(Complex, n_qreg)
        _n_qreg = n_qreg
        return Input(n_qreg, _init_state, [])
    end

    function apply!(circuit::QuantumCircuit, gate::Gate)
        push!(circuit._gates, gate)
        return circuit
    end

    function apply!(circuit::QuantumCircuit, gates::Array{T,1}) where T <: Gate
        for gate in gates
            push!(circuit._gates, gate)
        end
        return circuit
    end

    function apply!(circuit::QuantumCircuit, gates::Array{Any,1})
        for gate in gates
            if typeof(gate) <: Gate
                push!(circuit._gates, gate)
            else
                error("Incorrect gates.")
            end
        end
        return circuit
    end

    function apply!(circuit::Input, gate::Gate)
        push!(circuit._gates, gate)
        return circuit
    end

    function apply!(circuit::Input, gates::Array{T,1}) where T <: Gate
        for gate in gates
            push!(circuit._gates, gate)
        end
        return circuit
    end

    function apply!(circuit::Input, gates::Array{Any,1})
        for gate in gates
            if typeof(gate) <: Gate
                push!(circuit._gates, gate)
            else
                error("Incorrect gates.")
            end
        end
        return circuit
    end

    function execute(circuit::QuantumCircuit, device::T) where T <: Device
        result = execute_backend(circuit._n_qreg::Int, circuit._gates::Array{Gate,1}, device)
        return result
    end

    function execute(circuit::Input, device::T) where T <: Device
        result = execute_backend(circuit._n_qreg::Int, circuit._gates::Array{Gate,1}, device)
        return result
    end

end