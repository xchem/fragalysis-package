import argparse
import os

from frag.utils.network_utils import get_driver,add_attr,add_edge,add_node
from tqdm import tqdm


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Upload Nodes, Edges and Attributes')
    parser.add_argument('--base_dir')
    args = parser.parse_args()
    driver = get_driver()

    with driver.session() as session:
        for line in tqdm(open(os.path.join(args.base_dir,"nodes.txt")).readlines()):
            session.write_transaction(add_node,line.split()[1],line.split()[2],line.split()[3],line.split()[4])
        for line in tqdm(open(os.path.join(args.base_dir,"edges.txt")).readlines()):
            session.write_transaction(add_edge,line.split()[1],line.split()[2],line.split()[3])
        for line in tqdm(open(os.path.join(args.base_dir,"attributes.txt")).readlines()):
            session.write_transaction(add_attr,line.split()[1],line.split()[3])