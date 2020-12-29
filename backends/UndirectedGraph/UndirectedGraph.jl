module UndirectedGraphBackend

    include("./UndirectedGraphFunctions.jl")
    include("../backend_Models.jl")
    export execute_backend
    using .UndirectedGraphFunctions
    using .BackendModels

    function execute_backend(n_qregs::Int, gates, model)
        backend_gates = call_backend_gates(gates)
        parsed_gates = parse_gates(backend_gates)
        _gates = convert(Array{UG_Gate,1}, parsed_gates)
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

    function call_backend_gates(gates)
        backend_gates = []
        for gate in gates
            name = gate._name
            fields = fieldnames(typeof(gate))
            field_names = []
            for i in fields
                push!(field_names, string(i))
            end
            given_property = [getproperty(gate, property) for property in fields]
            code = gate._name * "("
            for prop in given_property[1:end-1]
                code *= string(prop)
                code *= ", "
            end
            code *= ")"
            new_gate = eval(Meta.parse(code))
            push!(backend_gates, new_gate)
        end
        return backend_gates
    end

    function parse_gates(backend_gates)
        parsed_gates =[]
        for gate in backend_gates
            _gates = gate_parser(gate)
            for _gate in _gates
                push!(parsed_gates, _gate)
            end
        end
        return parsed_gates
    end
    

end