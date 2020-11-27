from schema import Schema, Or, Optional

# Taskfile schemes
expandable_sequence = Schema({"from": float, "to": float, "by": float})
numeric_sequence = Schema({"start": float, "end": float, "step": float})

task_parameter = Or(float, int, str, expandable_sequence)

geometry_data = Schema({
    "type": str,
    "parameters": {str: task_parameter}
})

regular_metropolist_parameters = {
    Or("beta", "T", only_one=True): task_parameter,
    str: task_parameter
}

tempering_metropolis_parameters = {
    "tempering_betas": numeric_sequence,
    str: task_parameter
}

task = Schema({
    Optional("repeat"): int,
    "geometry": geometry_data,
    "parameters": Or(regular_metropolist_parameters, tempering_metropolis_parameters)
})

task_file_scheme = Schema({
    "pool_size": int,
    "output_folder": str,
    "tasks": [task]
})
