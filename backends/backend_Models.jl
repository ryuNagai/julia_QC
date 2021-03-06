module BackendModels

    export Device, Results
    export UndirectedGraphModel, StateVectorModel
    export StateVectorResults

    abstract type Device end
    abstract type Results end

    struct UndirectedGraphModel <: Device
        measure_all::Bool
        output_basis::Array{Int, 1}
    end

    function UndirectedGraphModel(x)
        if length(x) == 0
            return UndirectedGraphModel(true, [1])
        elseif typeof(x[1]) <: Array{Int, 1}
            for i in x[1]
                if i != 0 && i != 1
                    error("Measured state must be computational basis.")
                end
            end
            return UndirectedGraphModel(false, x[1])
        end
    end

    struct StateVectorModel <: Device
        zero_state::Bool
        init_state::Array{ComplexF64, 1}
    end

    function StateVectorModel(init_state)
        if length(init_state) == 0
            return StateVectorModel(true, [0im])
        else
            return StateVectorModel(false, init_state)
        end
    end

    struct StateVectorResults <: Results
        states::Array{Any, 1}
        cregs::Array{Any, 1}
    end
end