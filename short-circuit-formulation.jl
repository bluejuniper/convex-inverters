using PowerModelsProtection, PowerModelsDistribution, Ipopt, JuMP, Printf

const _PMD = PowerModelsDistribution
const _PMP = PowerModelsProtection

"Builds a multiconductor (distribution) fault study optimization problem"
function my_build_mc_fault_study(pm::_PMD.AbstractUnbalancedPowerModel)
    @debug "Building fault study"
    _PMD.variable_mc_bus_voltage(pm, bounded=false)
    _PMD.variable_mc_switch_current(pm, bounded=false)
    _PMD.variable_mc_branch_current(pm, bounded=false)
    _PMD.variable_mc_transformer_current(pm, bounded=false)
    _PMD.variable_mc_generator_current(pm, bounded=false)

    _PMP.variable_mc_bus_fault_current(pm)
    _PMP.variable_mc_pq_inverter(pm)
    _PMP.variable_mc_grid_formimg_inverter(pm)

    for (i,bus) in _PMD.ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        _PMD.constraint_mc_theta_ref(pm, i)
        _PMD.constraint_mc_voltage_magnitude_only(pm, i)
    end

    for id in _PMD.ids(pm, :gen)
        _PMD.constraint_mc_generator_power(pm, id; bounded=false)
    end

    # TODO add back in the generator voltage drop with inverters in model
    @debug "Adding constraints for synchronous generators"
    _PMP.constraint_mc_gen_voltage_drop(pm)

    for i in _PMD.ids(pm, :fault)
        _PMP.constraint_mc_bus_fault_current(pm, i)
    end

    for (i,bus) in _PMD.ref(pm, :bus)
        _PMP.constraint_mc_current_balance(pm, i)
    end

    for i in _PMD.ids(pm, :branch)
        _PMD.constraint_mc_current_from(pm, i)
        _PMD.constraint_mc_current_to(pm, i)
        _PMD.constraint_mc_bus_voltage_drop(pm, i)
        _PMP.expression_mc_branch_fault_sequence_current(pm, i)
    end

    for i in _PMD.ids(pm, :switch)
        _PMD.constraint_mc_switch_state(pm, i)
    end

    for i in _PMD.ids(pm, :transformer)
        _PMD.constraint_mc_transformer_power(pm, i)
    end

    @debug "Adding constraints for grid-following inverters"
    for i in _PMD.ids(pm, :solar_gfli)
        @debug "Adding constraints for grid-following inverter $i"
        _PMP.constraint_mc_pq_inverter(pm, i)
    end

    @debug "Adding constraints for grid-forming inverters"
    for i in _PMD.ids(pm, :solar_gfmi)
        @debug "Adding constraints for grid-forming inverter $i"
        # constraint_mc_grid_forming_inverter(pm, i)
        _PMP.constraint_mc_grid_forming_inverter_virtual_impedance(pm, i)
    end
end

function my_solve_mc_fault_study(case::Dict{String,<:Any}, solver; kwargs...)
    data = deepcopy(case)

    # TODO can this be moved?
    _PMP.check_microgrid!(data)

    solution = _PMD.solve_mc_model(
        data,
        _PMD.IVRUPowerModel,
        solver,
        build_mc_fault_study;
        eng2math_extensions=[_eng2math_fault!],
        eng2math_passthrough=_pmp_eng2math_passthrough,
        make_pu_extensions=[_PMP._rebase_pu_fault!, _PMP._rebase_pu_gen_dynamics!],
        map_math2eng_extensions=Dict{String,Function}("_map_math2eng_fault!"=>_map_math2eng_fault!),
        make_si_extensions=[_PMP.make_fault_si!],
        dimensionalize_math_extensions=_pmp_dimensionalize_math_extensions,
        ref_extensions=[_PMP.ref_add_mc_fault!, _PMP.ref_add_mc_solar!, _PMP.ref_add_grid_forming_bus!],
        solution_processors=[_PMP.solution_fs!],
        kwargs...
    )

    return solution
end

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