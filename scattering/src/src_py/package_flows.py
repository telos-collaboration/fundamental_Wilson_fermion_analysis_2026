#!/usr/bin/env python3

from argparse import ArgumentParser
import re

from flow_analysis.readers import read_flows_hirep
from flow_analysis.measurements.scales import measure_w0
import h5py
import numpy as np


def get_args():
    parser = ArgumentParser(
        description=(
            "Get gradient flow data from log files,"
            "and collate them into a single HDF5 file."
        )
    )

    parser.add_argument(
        "flow_filenames",
        nargs="+",
        metavar="flow_filename",
        help="Filename of gradient flow log file",
    )
    parser.add_argument(
        "--h5_filename",
        required=True,
        help="Where to place the combined HDF5 file.",
    )
    parser.add_argument(
        "--ensemble",
        required=True,
        help="Label of the ensmble to be used as the HDF5 groupname.",
    )
    parser.add_argument(
        "--W0",
        required=False,
        type=float,
        default=0.2815,
        help="Wilson flow reference scale.",
    )
    return parser.parse_args()


def parse_cfg_filename(filename):
    return re.match(
        r".*/(?P<runid>[^/]*)_(?P<NT>[0-9]+)x(?P<NX>[0-9]+)x(?P<NY>[0-9]+)x"
        r"(?P<NZ>[0-9]+)nc(?P<Nc>[0-9]+)(?:r(?P<rep>[A-Z]+))?(?:nf(?P<nf>[0-9]+))?"
        r"(?:b(?P<beta>[0-9]+\.[0-9]+))?(?:m(?P<mass>-?[0-9]+\.[0-9]+))?"
        r"n(?P<cfg_idx>[0-9]+)",
        filename,
    ).groupdict()


def get_filename_metadata(metadata, content):
    if content[0] == "[IO][0]Configuration" and content[2] == "read":
        provisional_metadata = parse_cfg_filename(content[1])
        keys = {
            "NT": int,
            "NX": int,
            "NY": int,
            "NZ": int,
            "Nc": int,
            "mass": float,
            "beta": float,
        }

        for key, dtype in keys.items():
            value = dtype(provisional_metadata[key])
            if key in metadata and metadata[key] != value:
                message = f"Metadata inconsistent: {metadata[key]} != {value}."
                raise ValueError(message)
            metadata[key] = value


def process_file(flow_filename, h5file, W0, group_name):
    flows = read_flows_hirep(flow_filename, metadata_callback=get_filename_metadata)
    group = h5file.create_group(group_name)
    group.create_dataset("beta", data=flows.metadata["beta"])
    group.create_dataset("configurations", data=flows.cfg_filenames.astype("S"))
    group.create_dataset("trajectory indices", data=flows.trajectories)
    group.create_dataset("ensemble names", data=flows.ensemble_names.astype("S"))
    group.create_dataset("gauge group", data=f"SP({flows.metadata['Nc']})")
    group.create_dataset(
        "lattice",
        data=np.asarray([flows.metadata[key] for key in ["NT", "NX", "NY", "NZ"]]),
    )
    group.create_dataset("plaquette", data=flows.plaquettes)
    group.create_dataset("quarkmass", data=[flows.metadata["mass"]])
    group.create_dataset("flow type", data=flows.metadata.get("flow_type"))

    group.create_dataset("flow times", data=flows.times)
    group.create_dataset("topological charge", data=flows.Qs)
    group.create_dataset("energy density plaq", data=flows.Eps)
    group.create_dataset("energy density sym", data=flows.Ecs)

    # add additional quantities for compatibility
    Qs = flows.Q_history()
    trajectories = flows.trajectories
    w0_sym = measure_w0(flows, W0, operator="sym")
    w0_plaq = measure_w0(flows, W0, operator="plaq")

    # (I follow the analysis release of )
    flow_time_index_sym  = abs(flows.times - w0_sym**2).argmin()
    flow_time_index_plaq = abs(flows.times - w0_plaq**2).argmin()
    energy_density_w0_sym = flows.Ecs[:,flow_time_index_sym]
    energy_density_w0_plaq = flows.Eps[:,flow_time_index_plaq]
    group.create_dataset("trajectories", data=trajectories)
    group.create_dataset("Q", data=Qs)
    group.create_dataset("w0_val", data=w0_sym.nominal_value)
    group.create_dataset("w0_std", data=w0_sym.std_dev)
    group.create_dataset("energy_density_w0_sym",  data=energy_density_w0_sym)
    group.create_dataset("energy_density_w0_plaq", data=energy_density_w0_plaq)


def main():
    args = get_args()
    with h5py.File(args.h5_filename, "a") as h5file:
        for flow_filename in args.flow_filenames:
            process_file(flow_filename, h5file, args.W0, args.ensemble)

if __name__ == "__main__":
    main()
