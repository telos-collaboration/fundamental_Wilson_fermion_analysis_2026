for filename in raw_data/raw_flow/*/topology/out/out_flow; do 
    name=$(echo $filename | grep -o "Lt.*FUN")
    echo $name 
    python3 src/src_py/package_flows.py --h5_filename data_assets/topology.hdf5 --ensemble $name $filename
done