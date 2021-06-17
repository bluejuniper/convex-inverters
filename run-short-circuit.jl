using PowerModelsProtection, Ipopt, JuMP, Printf

net = PowerModelsProtection.parse_file("case3_unbalanced.dss")
net["multinetwork"] = false
solver = JuMP.with_optimizer(Ipopt.Optimizer)

# Simulate the fault
net["fault"] = Dict{String, Any}()
PowerModelsProtection.add_fault!(net, "1", "lg", "loadbus", [1, 4], 0.005)
results = PowerModelsProtection.solve_mc_fault_study(net, solver)

# Print out the fault currents
Iabc = results["solution"]["line"]["ohline"]["fault_current"]
@printf("Fault current: %0.3f A\n", Iabc[1])