# import packages
import os
import h5py
import shortuuid
from schema import Schema, Or
import json
from itertools import product
from functools import reduce
import numpy as np

# local scripts
from scripts.mpi_logger import get_logger
from . import schemes

# Get MPI communicator
logger = get_logger('io')
logger.debug('Starting io module')

metadata_scheme = Schema({
    "geometry": {"type": str, "parameters": dict},
    "parameters": dict,
    "output_folder": str
})

numerical_filter = Schema(Or(int, float))


class SpintronicsSnapshoter:
    def __init__(self, spintronics, filename, metadata):
        logger.debug(f"Starting hdf5 snapshoter for file: {filename}")

        self.spintronics = spintronics
        self.file = h5py.File(filename, 'a')
        self.filename = filename

        logger.debug(f"Writing metadata as attributes")
        try:
            metadata_scheme.validate(metadata)
        except Exception as e:
            logger.error(f"Invalid attributes scheme: {getattr(e, 'message', str(e))}")

        numerical_filter = Schema(Or(int, float))
        numerical_geometry_data = {k: v for (
            k, v) in metadata['geometry']['parameters'].items() if numerical_filter.is_valid(v)}

        # Set attributes
        self.file.attrs['geometry-type'] = metadata['geometry']['type']
        self.file.attrs['geometry-attrs'] = list(
            numerical_geometry_data.keys())
        self.file.attrs['geometry-parms'] = list(
            numerical_geometry_data.values())

    def snapshot(self, **kwargs):
        n = len(self.spintronics.rx)
        snapshot_name = f"snap-{shortuuid.uuid()}"

        #  check colisions (just to be sure)
        while snapshot_name in list(self.file.keys()):
            snapshot_name = f"snap-{shortuuid.uuid()}"

        snapshot = self.file.create_dataset(snapshot_name, (n, 3), dtype='f')
        snapshot[:, 0] = self.spintronics.rx
        snapshot[:, 1] = self.spintronics.ry
        snapshot[:, 2] = self.spintronics.rz

        # logger.debug(f"Took snapshot with size {snapshot.shape}")

        for (k, v) in kwargs.items():
            snapshot.attrs[k] = v

    def update_attrs(self, **kwargs):
        for (k, v) in kwargs.items():
            self.file.attrs[k] = v

    def close(self):
        self.file.close()
        logger.debug(f"Closing file: {self.filename}")


class TaskIO:
    def __init__(self, filefolder, filename):
        logger.debug("Initializing resume writter")
        self.filename = filename
        self.filefolder = filefolder
        self.current_dset = None

    @staticmethod
    def save_record(file, record):
        record_geometry = record['geometry']
        record_geometry['parameters'] = {key: value for (
            key, value) in record_geometry['parameters'].items() if numerical_filter.is_valid(value)}

        geometry_dset = None
        for dset_name in file.keys():
            dset_geometry_parms_zip = zip(
                file[dset_name].attrs['geometry-attrs'], file[dset_name].attrs['geometry-parms'])

            dset_geometry_data = {
                "type": file[dset_name].attrs['geometry-type'],
                "parameters": {x[0]: x[1] for x in dset_geometry_parms_zip}
            }

            if(dset_geometry_data == record_geometry):
                geometry_dset = file[dset_name]
                break

        record_data = {key: value for (
            key, value) in record.items()
            if key != 'geometry' and numerical_filter.is_valid(value)}

        if(geometry_dset == None):
            logger.debug(f"Creating dataset for geometry {record['geometry']}")

            geometry_dset_name = f"resumes-{shortuuid.uuid()}"
            n_cols = len(record_data)

            geometry_dset = file.create_dataset(
                geometry_dset_name, (0, n_cols), dtype='f', maxshape=(None, n_cols))

            # Write metadata
            geometry_dset.attrs['columns'] = sorted(record_data.keys())

            numerical_geometry_data = {k: v for (
                k, v) in record['geometry']['parameters'].items() if numerical_filter.is_valid(v)}

            geometry_dset.attrs['geometry-attrs'] = list(
                numerical_geometry_data.keys())
            geometry_dset.attrs['geometry-parms'] = list(
                numerical_geometry_data.values())
            geometry_dset.attrs['geometry-type'] = record['geometry']['type']

        # Reshape for one more entry
        columns = geometry_dset.attrs['columns']
        n_rows = geometry_dset.shape[0]
        geometry_dset.resize(n_rows + 1, axis=0)
        row_data = [record_data[column] for column in columns]

        geometry_dset[n_rows, :] = row_data

    def save(self, records):
        logger.debug(f"Writing {len(records)} to file {self.filename}")

        filepath = os.path.join(self.filefolder, self.filename)
        with h5py.File(filepath, 'a') as f:
            for record in records:
                TaskIO.save_record(f, record)

        logger.debug(f"Writing {len(records)} to file ~{self.filename}")

        # Save a safe copy
        filepath = os.path.join(self.filefolder, f"~{self.filename}")
        with h5py.File(filepath, 'a') as f:
            for record in records:
                TaskIO.save_record(f, record)


def expand_parameter_in_sequence(parameter):
    if schemes.expandable_sequence.is_valid(parameter) or schemes.numeric_sequence.is_valid(parameter):
        n = max(0, np.rint(
            (parameter['to'] - parameter['from'])/parameter['by']))
        n = int(n)

        if(parameter['from'] + n*parameter['by'] <= parameter['to']):
            n += 1

        if(schemes.numeric_sequence.is_valid(parameter)):
            return [np.linspace(parameter['from'], parameter['to'], n)]
        else:
            return np.linspace(parameter['from'], parameter['to'], n)
    else:
        return [parameter]


def merge(acc, cur):
    return dict(acc, **cur)


def read_input_tasks(filepath):
    data = None
    with open(filepath, 'r') as f:
        data = json.load(f)

    try:
        logger.info(f"Validating task file: {filepath}")
        schemes.task_file_scheme.validate(data)
        logger.info("Task file is valid")
    except Exception as e:
        logger.error(f"Invalid task file: {filepath}")
        logger.error(f"Scheme error {getattr(e, 'message', str(e))}")
        raise e

    tasks = []

    for task_data in data['tasks']:
        geometry_parms = [
            [{k: v} for v in expand_parameter_in_sequence(parm)] for k, parm in task_data['geometry']["parameters"].items()
        ]

        task_parms = [
            [{k: v} for v in expand_parameter_in_sequence(parm)] for k, parm in task_data["parameters"].items()
        ]

        repeat = task_data['repeat'] if 'repeat' in task_data else 1
        for geom_parm_point in product(*geometry_parms):
            for task_parms_point in product(*task_parms):

                task = {
                    "output_folder": data['output_folder'],
                    "geometry": {
                        "type": task_data['geometry']["type"],
                        "parameters": reduce(merge, [task_data['geometry']["parameters"], *geom_parm_point])
                    },
                    "parameters": reduce(merge, [task_data["parameters"], *task_parms_point])
                }
                for r in range(repeat):
                    tasks.append(task)
    return {"tasks": tasks, "pool_size": data["pool_size"]}


def get_task_filepath(output_folder):
    task_data_filename = f"task-{shortuuid.uuid()}.h5"
    task_data_filepath = os.path.join(output_folder, task_data_filename)

    # Ensure colision
    while(os.path.exists(task_data_filepath)):
        task_data_filename = f"task-{shortuuid.uuid()}.h5"
        task_data_filepath = os.path.join(output_folder, task_data_filename)

    return task_data_filepath
