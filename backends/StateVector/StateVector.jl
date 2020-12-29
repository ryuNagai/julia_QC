module StateVectorBackend

    include("./StateVectorFunctions.jl")
    include("../backend_Models.jl")
    export execute_backend
    using .StateVectorFunctions
    using .BackendModels

    function execute_backend(n_qregs::Int, gates, model)
        if model.zero_state == true
            state = zeros(Complex, 2^n_qregs)
            state[1] = 1.
        else
            if length(model.init_state) != 2^n_qregs
                error("Initial state length must be equal to 2^(n_qubits).")
            else
                state = model.init_state
            end
        end
        backend_gates = call_backend_gates(gates)
        for gate in backend_gates
            state = apply!(gate, state, n_qregs)
        end
        return state
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