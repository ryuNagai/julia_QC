module StateVectorBackend

    include("./StateVectorFunctions.jl")
    include("../backend_Models.jl")
    export execute_backend
    using .StateVectorFunctions
    using .BackendModels

    function execute_backend(n_qregs::Int, gates, model, counts::Int64)
        if model.zero_state == true
            state = zeros(ComplexF64, 2^n_qregs)
            state[1] = 1.
        else
            if length(model.init_state) != 2^n_qregs
                error("Initial state length must be equal to 2^(n_qubits).")
            else
                state = model.init_state
            end
        end

        backend_gates = call_backend_gates(gates)
        states = []
        cregs = []
        for i in 1:counts
            _backend_gates = copy(backend_gates)
            _state = copy(state)
            creg_indices = []
            results = []
            for gate in _backend_gates
                apply!(gate, _state, n_qregs)
                if typeof(gate) == M
                    push!(creg_indices, gate._cregidx)
                    push!(results, gate._result)
                end
            end

            if length(creg_indices) > 0
                creg = zeros(Int64, max(creg_indices...))
                for (i, idx) in enumerate(creg_indices)
                    creg[idx] = results[i]
                end
            else
                creg = zeros(Int64, n_qregs)
            end
            push!(states, _state)
            push!(cregs, creg)
        end
        return states, cregs
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

end