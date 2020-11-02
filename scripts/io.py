# import packages
import os
import h5py
import shortuuid
from schema import Schema, Or
import json

# local scripts
from scripts.mpi_logger import get_logger

# Get MPI communicator
logger = get_logger('io')
logger.debug('Starting io module')

metadata_scheme = Schema({
    "geometry": {"type": str, "parameters": dict},
    "parameters": dict
})


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

        logger.debug(f"Took snapshot with size {snapshot.shape}")

        for (k, v) in kwargs.items():
            snapshot.attrs[k] = v

    def close(self):
        self.file.close()
        logger.debug(f"Closing file: {self.filename}")


class TaskIO:
    def __init__(self, filepath):
        logger.debug("Initializing resume writter")
        self.filepath = filepath
        self.current_dset = None

    @staticmethod
    def save_record(file, record):
        geometry_dset = None
        for dset_name in file.keys():
            dset_geometry_parms_zip = zip(
                file[dset_name].attrs['geometry-attrs'], file[dset_name].attrs['geometry-parms'])

            dset_geometry_data = {
                "type": file[dset_name].attrs['geometry-type'],
                "parameters": {x[0]: x[1] for x in dset_geometry_parms_zip}
            }

            if(dset_geometry_data == record["geometry"]):
                geometry_dset = file[dset_name]
                break

        record_data = {key: value for (
            key, value) in record.items() if key != 'geometry'}

        if(geometry_dset == None):
            geometry_dset_name = f"resumes-{shortuuid.uuid()}"
            n_cols = len(record_data)

            geometry_dset = file.create_dataset(
                geometry_dset_name, (0, n_cols), dtype='f', maxshape=(None, n_cols))

            # Write metadata
            geometry_dset.attrs['columns'] = sorted(record_data.keys())

            geometry_dset.attrs['geometry-attrs'] = list(
                record['geometry']['parameters'].keys())
            geometry_dset.attrs['geometry-parms'] = list(
                record['geometry']['parameters'].values())
            geometry_dset.attrs['geometry-type'] = record['geometry']['type']

        # Reshape for one more entry
        columns = geometry_dset.attrs['columns']
        n_rows = geometry_dset.shape[0]
        geometry_dset.resize(n_rows + 1, axis=0)
        row_data = [record_data[column] for column in columns]

        geometry_dset[n_rows, :] = row_data

    def save(self, records):
        logger.debug(f"Writing {len(records)} to file {self.filepath}")

        with h5py.File(self.filepath, 'a') as f:
            for record in records:
                TaskIO.save_record(f, record)

        logger.debug(f"Writing {len(records)} to file ~{self.filepath}")

        # Save a safe copy
        with h5py.File(f"~{self.filepath}", 'a') as f:
            for record in records:
                TaskIO.save_record(f, record)


def read_input_tasks(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)

        return data
