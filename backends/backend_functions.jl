module BackendFunctions
    include("./UndirectedGraph/UndirectedGraph.jl")
    include("./StateVector/StateVector.jl")
    include("./backend_Models.jl")
    include("../Gate.jl")

    import .UndirectedGraphBackend
    import .StateVectorBackend
    using Reexport
    @reexport using .GateSet
    @reexport using .BackendModels

    export execute_backend

    function execute_backend(n_qregs::Int, gates::Array{Gate,1}, model::UndirectedGraphModel)
        return UndirectedGraphBackend.execute_backend(n_qregs::Int, gates, model)
    end

    function execute_backend(n_qregs::Int, gates::Array{Gate,1}, model::StateVectorModel)
        return StateVectorBackend.execute_backend(n_qregs::Int, gates, model)
    end
end