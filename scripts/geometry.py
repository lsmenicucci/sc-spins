# import packages
import logging
import numpy as np
import json
from mpi4py import MPI
from schema import Schema, And

# import local modules
from scripts.mpi_logger import get_logger

# Get MPI communicator
logger = get_logger('geometry')
logger.debug('Starting geometry module')

parameters_schemas = {
    "disk": Schema({
        "r": int,
        "layers": int,
        "a": float,
        "dipolar_cut": float,
        "J": float,
        "D": float
    }),
    "square": Schema({
        "L": int,
        "layers": int,
        "a": float,
        "dipolar_cut": float,
        "J": float,
        "D": float
    }),
    "line": Schema({
        "L": int,
        "a": float,
        "dipolar_cut": float,
        "J": float,
        "D": float
    }),
    "kekulene": Schema({
        "filename": str,
        "a": float,
        "dipolar_cut": float,
        "J": float,
        "D": float,
        "S_carbon": float,
        "S_hydrogen": float
    })
}


def get_sites_in_range(spintronics, center_site, cut_range, include_border=False):
    n = len(spintronics.rx)

    site_x = spintronics.rx[center_site]
    site_y = spintronics.ry[center_site]
    site_z = spintronics.rz[center_site]

    def dist(i):
        return np.sqrt(
            (site_x - spintronics.rx[i])**2 +
            (site_y - spintronics.ry[i])**2 +
            (site_z - spintronics.rz[i])**2)

    if(include_border):
        return list(filter(lambda s: s != center_site and dist(s) <= cut_range, range(n)))
    else:
        return list(filter(lambda s: s != center_site and dist(s) < cut_range, range(n)))


def initialize_general_geometry(spintronics, points, a=1.0, dipolar_cut=None, J=1.0, D=1.0, **kwargs):
    logger.debug(f'Initializin general geometry with {len(points)} points ')
    n = len(points)

    # Set the sites positions
    spintronics.rx = list(map(lambda p: p[0], points))
    spintronics.ry = list(map(lambda p: p[1], points))
    spintronics.rz = list(map(lambda p: p[2], points))

    # Initialize spin vectors
    spintronics.sx = np.zeros(n)
    spintronics.sy = np.zeros(n)
    spintronics.sz = np.ones(n)

    # Initialize spin constants
    spintronics.spin = np.ones(n)

    if('spin' in kwargs):
        spintronics.spin = kwargs["spin"]

    # get neibs
    if('J_neibs' in kwargs):
        logging.debug(f"Using provided exchange neibs with size = {len(kwargs['J_neibs'])}")
        exc_neibs = kwargs["J_neibs"]
    else:
        exc_neibs = list(map(lambda s: get_sites_in_range(
            spintronics, s, a, include_border=True), range(n)))
    dip_neibs = list(map(lambda s: get_sites_in_range(
        spintronics, s, dipolar_cut), range(n)))

    max_exc_neibs_count = max(map(lambda neibs: len(neibs), exc_neibs))
    max_exc_dip_count = max(map(lambda neibs: len(neibs), dip_neibs))

    # Initalize exchange interaction vectors
    spintronics.v_exc = J * np.ones([3, max_exc_neibs_count, n])
    spintronics.v_interacts_exc = -1 * np.ones([max_exc_neibs_count, n])
    spintronics.v_interacts_exc_count = -1 * np.ones(n)

    # Initialize dipolar interaction vectors
    spintronics.v_dip = D * np.ones([3, max_exc_dip_count, n])
    spintronics.v_interacts_dip = -1 * np.ones([max_exc_dip_count, n])
    spintronics.v_interacts_dip_count = -1 * np.ones(n)

    logger.debug('Calculating neibs')

    # Load neib data
    for i in range(n):
        exc_neib_size = len(exc_neibs[i])
        dip_neib_size = len(dip_neibs[i])

        spintronics.v_interacts_exc_count[i] = exc_neib_size
        for exc_j in range(exc_neib_size):
            spintronics.v_interacts_exc[exc_j, i] = exc_neibs[i][exc_j] + 1

        spintronics.v_interacts_dip_count[i] = dip_neib_size
        for dip_j in range(dip_neib_size):
            spintronics.v_interacts_dip[dip_j, i] = dip_neibs[i][dip_j] + 1

    logger.debug('Done setting geometry')
    return points


def initialize_disk(spintronics, parameters):
    try:
        parameters_schemas["disk"].validate(parameters)
    except Exception as e:
        logger.error(f"Error validating geometry parameters: {getattr(e, 'message', str(e))}")
        raise

    logger.info(
        f'Initializing disk geometry with {parameters["layers"]} layers and radius = {parameters["r"]}')

    slimit = parameters["r"] * int(np.round(float(parameters["r"]) / 2.0))
    points = [(i * parameters["a"], j * parameters["a"], k * parameters["a"])
              for k in range(parameters["layers"])
              for j in range(-slimit, slimit + 1)
              for i in range(-slimit, slimit + 1)
              if i**2 + j**2 < parameters["r"]**2]

    initialize_general_geometry(spintronics, points, **parameters)
    return {"type": "disk", "parameters": parameters}


def initialize_square_layers(spintronics, parameters):
    try:
        parameters_schemas["square"].validate(parameters)
    except Exception as e:
        logger.error(f"Error validating geometry parameters: {getattr(e, 'message', str(e))}")
        raise

    logger.info(
        f'Initializing square geometry with {parameters["layers"]} layers and side = {parameters["L"]}')

    slimit = int(np.round(float(parameters["L"]) / 2.0))
    points = [(i * parameters["a"], j * parameters["a"], k * parameters["a"])
              for k in range(parameters["layers"])
              for j in range(-slimit, slimit + 1)
              for i in range(-slimit, slimit + 1)]

    initialize_general_geometry(spintronics, points, **parameters)
    return {"type": "square", "parameters": parameters}


def initialize_line(spintronics, parameters):
    try:
        parameters_schemas["line"].validate(parameters)
    except Exception as e:
        logger.error(f"Error validating geometry parameters: {getattr(e, 'message', str(e))}")
        raise

    logger.info(f'Initializing line geometry with length = {parameters["L"]}')

    points = [(i * parameters["a"], 0, 0)
              for i in range(-int(parameters["L"] / 2), int(parameters["L"] / 2))]

    initialize_general_geometry(spintronics, points, **parameters)
    return {"type": "line", "parameters": parameters}


def initialize_organic_molecule(spintronics, parameters):
    try:
        parameters_schemas["kekulene"].validate(parameters)
    except Exception as e:
        logger.error(f"Error validating geometry parameters: {getattr(e, 'message', str(e))}")
        raise

    logger.info(f"Initializing organic molecule geometry")
    logger.info(f"Opening data file: {parameters['filename']}")

    data = None
    with open(parameters["filename"], 'r') as datafile:
        data = json.load(datafile)

    # Per atom-data
    points = []
    J_neibs = []
    spin = []

    n_carbons = len(data['carbon_atoms'])
    for (n, atom) in data['carbon_atoms'].items():
        neibs = [(neib - 1) + (neib_typ * n_carbons)
                 for (neib, neib_typ) in zip(atom['neighbors'], atom['type'])]

        points.append(atom['position'] + [0.0])
        J_neibs.append(neibs)
        spin.append(parameters["S_carbon"])

    for (n, atom) in data['hydrogen_atoms'].items():
        neibs = [(neib - 1) for neib in atom['neighbors']]

        points.append(atom['position'] + [0.0])
        J_neibs.append(neibs)
        spin.append(parameters["S_hydrogen"])

    initialize_general_geometry(
        spintronics, points, J_neibs=J_neibs, spin=spin, **parameters)

    return {"type": "kekulene", "parameters": parameters}


def initialize(spintronics, geometryType, parameters):
    logger.info("Initializing geometry")

    if(geometryType == 'disk'):
        return initialize_disk(spintronics, parameters)
    elif(geometryType == 'square'):
        return initialize_square_layers(spintronics, parameters)
    elif(geometryType == 'line'):
        return initialize_line(spintronics, parameters)
    elif(geometryType == 'kekulene'):
        return initialize_organic_molecule(spintronics, parameters)
    else:
        logger.error(f"Can't initialize geometry named: ${geometryType}")
        raise RuntimeError(f'Undefined geometry type: {geometryType}')
