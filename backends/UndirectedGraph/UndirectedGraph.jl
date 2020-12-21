module UndirectedGraphBackend

    include("./UndirectedGraphFunctions.jl")
    export UndirectedGraphModel, execute_backend
    using Reexport
    @reexport using .UndirectedGraphFunctions

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


    function execute_backend(n_qregs::Int, gates::Array{Gate,1}, model::UndirectedGraphModel)
        _gates =[]
        for gate in gates
            parsed_gates = gate_parser(gate)
            _gates = vcat(_gates, parsed_gates)
        end
        if model.measure_all
            input_state = zeros(Int, n_qregs)
            output_state = zeros(Int, n_qregs)
            state_vec = []
            buckets = create_buckets(n_qregs, _gates, input_state, output_state)
            result = contract_graph(buckets)
            push!(state_vec, result)
            while true
                for i in 1:n_qregs
                    if output_state[n_qregs - (i-1)] == 0
                        output_state[n_qregs - (i-1)] = 1
                        output_state[n_qregs - (i-1) + 1:n_qregs] .= 0
                        break
                    end
                end
                buckets = create_buckets(n_qregs, _gates, input_state, output_state)
                result = contract_graph(buckets)
                push!(state_vec, result)
                if sum(output_state) >= n_qregs
                    break
                end
            end
            return state_vec
        else
            if length(model.output_basis) != n_qregs
                error("Measured state length must be equal to qubit number.")
            end
            input_state = zeros(Int, n_qregs)
            output_state = model.output_basis
            buckets = create_buckets(n_qregs, _gates, input_state, output_state)
            result = contract_graph(buckets)
            return result
        end
        
    end
    

end