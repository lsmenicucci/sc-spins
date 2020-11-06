from schema import Schema, Or, Optional

# Taskfile schemes
numeric_sequence = Schema({"from": float, "to": float, "by": float})

task_parameter = Or(float, int, str, numeric_sequence)

geometry_data = Schema({
    "type": str,
    "parameters": {str: task_parameter}
})

task = Schema({
    Optional("repeat"): int,
    "geometry": geometry_data,
    "parameters": {
        Or("beta", "T", only_one=True): task_parameter,
        str: task_parameter
    }
})

task_file_scheme = Schema({
    "pool_size": int,
    "output_folder": str,
    "tasks": [task]
})
