#!/usr/bin/env python3
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   kmtricks
#   Authors: T. Lemane
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys
import os
import argparse
import subprocess
from multiprocessing import cpu_count
import shutil
import pickle
from typing import Tuple, List, Dict, Any, NamedTuple
from kmtricks import __version__

class OutputFormatter:
    def __init__(self, content: str, map_id: dict, mode: str, path: str, query: str):
        self.content = content
        self.path = path
        self.map = map_id
        self.query_in = iter(open(query, 'r'))
        if mode == 'list':
            self.format = self.as_list
        else:
            self.format = self.as_bit_vector

    def as_bit_vector(self):
        with open(self.path, 'w') as out:
            name = ''
            res = ['0']*(len(self.map)//2)
            for line in self.content.split('\n'):
                if line:
                    if line.startswith('*'):
                        if name:
                            out.write(f'{name}: {"".join(res)}\n')
                            res = ['0']*(len(self.map)//2)
                            name = next(self.query_in).rstrip()
                        else:
                            name = next(self.query_in).rstrip()
                    else:
                        res[self.map[line]] = '1'
            out.write(f'{name}: {"".join(res)}\n')
    
    def as_list(self):
        with open(self.path, 'w') as out:
            name = ''
            res = []
            for line in self.content.split('\n'):
                if line:
                    if line.startswith('*'):
                        if name:
                            out.write(f'{name}: {" ".join(res)}\n')
                            res = []
                            name = next(self.query_in).rstrip()
                        else:
                            name = next(self.query_in).rstrip()
                    else:
                        res.append(line)
            out.write(f'{name}: {"".join(res)}\n')
                

class KmtricksCmd:
    cmd = (
      "kmtricks.py run --run-dir {output} --file {fof} --kmer-size {kmer_size} "
      "--count-abundance-min {min_abundance} --nb-cores {nb_cores} --hasher sabuhash --mode bf_trp "
      "--max-hash {bf_size} --split howde"
    )
    def __init__(self, args: argparse.Namespace):
        self.args = args
    
    def exec(self):
        fmt_cmd = KmtricksCmd.cmd.format(**vars(self.args))
        p = subprocess.Popen(fmt_cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
        if p.returncode != 0:
            sys.exit('kmtricks fails.')

class HowDeBuildCmd:
    cmd0 = "ls *.bf > index"
    cmd1 = "km_howdesbt cluster index"
    cmd2 = "km_howdesbt build --howde index"
    def __init__(self):
        pass

    def exec(self):
        subprocess.call(HowDeBuildCmd.cmd0, shell=True)
        p1 = subprocess.Popen(
          HowDeBuildCmd.cmd1.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p1.wait()
        if p1.returncode != 0:
            sys.exit('km_howdesbt cluster fails.')
        
        p2 = subprocess.Popen(
          HowDeBuildCmd.cmd2.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p2.wait()
        if p2.returncode != 0:
            sys.exit('km_howdesbt build fails.')

class HowDeQueryCmd:
    cmd = "km_howdesbt queryKm --tree={tree} --repart={repart} --win={win} {query}"
    def __init__(self, args: argparse.Namespace):
        self.args = args
    
    def exec(self):
        fmt_cmd = HowDeQueryCmd.cmd.format(**vars(self.args))
        p = subprocess.Popen(
          fmt_cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p.wait()
        if p.returncode != 0:
            sys.exit('km_howdesbt queryKm fails.')
        
        out, err = p.communicate()
        formatter = OutputFormatter(out.decode('utf8'),
                                    self.args.map,
                                    self.args.output_type,
                                    self.args.output,
                                    self.args.query)
        formatter.format()

def convert_to_kmtricks_fof(socks_input: str, fof_output: str) -> None:
    map_id = {}
    with open(socks_input, 'r') as f_in:
        with open(fof_output, 'w') as f_out:
            i = 0
            for line in f_in:
                line = line.rstrip()
                if line:
                    idx, paths = line.split(':')
                    map_id[i] = idx; map_id[idx] = i
                    i += 1
                    km_format = f'{idx}: {";".join(paths.strip().split(" "))}\n'
                    f_out.write(km_format)
    return map_id

class Cd: 
    def __init__(self, dirname):
        self.dirname = dirname

    def __enter__(self):
        self.current = os.getcwd()
        os.chdir(self.dirname)

    def __exit__(self, type, value, traceback):
        os.chdir(self.current)

SUBCOMMANDS = ('build', 'lookup-kmer')

def parser(description: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description)

    p._positionals.title = 'Subcommands'
    p.add_argument('--version', action='version',
      version=f'kmtricks-socks v{__version__}', help='Display kmtricks-socks version')
    subparsers = p.add_subparsers(
      dest='cmd', metavar='kmtricks-socks.py', help=f'[{" | ".join(SUBCOMMANDS)}]')

    build_parser = subparsers.add_parser('build', description="")
    build_parser.add_argument('-k', '--kmer-size', type=int, metavar='INT', default=20,
      help='Size of k-mers (default: %(default)s)')
    build_parser.add_argument('-m', '--min-abundance', type=int, metavar='INT', default=2,
      help='Min abundance for solid k-mers (default: %(default)s)')
    build_parser.add_argument('-b', '--bf-size', type=int, metavar='INT', default=10000,
      help='Bloom filter size (default: %(default)s)')
    build_parser.add_argument('-n', '--nb-cores', type=int, metavar='INT',
      default=cpu_count(), help='Number of cores (default: %(default)s)')
    build_parser.add_argument(dest='input', nargs=1, type=str)
    build_parser.add_argument(dest='output', nargs='?',
      default='./km-index', help='(default: %(default)s)')
    
    lookup_parser = subparsers.add_parser('lookup-kmer', description="")
    lookup_parser.add_argument(
      '-i', '--index-dir', type=str, metavar='STR', default='./km-index',
      help='(default: %(default)s)')
    lookup_parser.add_argument(
      '-t', '--threshold', type=float, metavar='FLOAT', default=1.0,
      help='(default: %(default)s)')
    lookup_parser.add_argument('-o', '--output-type', metavar='STR', type=str,
            choices=('vector', 'list'), default='vector',
            help=f'output type: [{"|".join(["vector", "list"])}] (default: %(default)s)')

    lookup_parser.add_argument('input', nargs=1)
    lookup_parser.add_argument(
      'output', nargs='?', default='results.txt', help='(default: %(default)s)')
    return p

def parse_args(parser) -> argparse.Namespace:
    if len(sys.argv) < 2:
        return parser.parse_args(['-h'])
    return parser.parse_args()
    

def main():
    description = 'kmtricks socks interface, see https://gitlab.ub.uni-bielefeld.de/gi/socks'
    args = parse_args(parser(description))

    if not shutil.which("kmtricks.py"):
        print('kmtricks.py not in path.')
        sys.exit(1)

    if args.cmd == 'build':
        km_fof = f'{"".join(args.input)}.kmtricks'
        map_id = convert_to_kmtricks_fof(''.join(args.input), km_fof)
        args.fof = km_fof
        km_cmd = KmtricksCmd(args)
        km_cmd.exec()

        with Cd(f'{args.output}/storage/vectors/howde/'):
            hb_cmd = HowDeBuildCmd()
            hb_cmd.exec()
        with open(f'{args.output}/socks.pickle', 'wb') as p_out:
            pickle.dump(map_id, p_out)
        print(f'Done. Query: kmtricks-socks.py lookup-kmer <query_file> -i {args.output}')

    elif args.cmd == 'lookup-kmer':
        if not os.path.exists(f'{args.index_dir}'):
            print("Index not found.")
            sys.exit(1)
        args.input = ''.join(args.input)
        args.tree = f'{args.index_dir}/storage/vectors/howde/index.detbrief.sbt'
        args.repart = f'{args.index_dir}/storage/partition_storage_gatb/minimRepart.minimRepart'
        args.win = f'{args.index_dir}/storage/hash_window.vec'
        args.query = args.input
        with open(f'{args.index_dir}/socks.pickle', 'rb') as p_in:
            args.map = pickle.load(p_in)
        args.mode = 2
        hq_cmd = HowDeQueryCmd(args)
        hq_cmd.exec()

        print(f'Results at {args.output}.')

if __name__ == '__main__':
    main()