module StateVectorBackend

    include("./StateVectorFunctions.jl")
    export StateVectorModel, execute_backend
    using Reexport
    @reexport using .StateVectorFunctions

    struct StateVectorModel <: Device
        init_state::Array{Complex, 1}
    end

    function StateVectorModel(init_state::Array{Complex, 1})
        return StateVectorModel(init_state)
    end

    function execute_backend(n_qregs::Int, gates::Array{Gate,1}, model::StateVectorModel)
        # _gates =[]
        # for gate in gates
        #     parsed_gates = gate_parser(gate)
        #     for _gate in parsed_gates
        #         push!(_gates, _gate)
        #     end
        # end
        # _gates = convert(Array{Gate,1}, _gates)
        state = model.init_state
        for gate in gates
            state = apply!(gate, state, n_qregs)
        end
        return state
    end
    

end