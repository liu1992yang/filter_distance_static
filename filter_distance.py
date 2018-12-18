import numpy as np
import pandas as pd
import math
import re
import sys
import os
from itertools import *
from scipy.spatial import distance

FILTER = 4 #Angstrong
wdir = '.'


def get_coords(atom_xyz,index):
    return atom_xyz[index][1:]

def dist(atom_xyz, index1, index2):
  return distance.euclidean(get_coords(atom_xyz,index1), get_coords(atom_xyz,index2))

def read_onefile(filename):
  with open(os.path.join(wdir, filename),'r') as fin:
    data=[]
    for line in fin:
      if line.strip() =='':
        continue
      if line.startswith("  "):
        arr = re.split(" +", line.strip())
        if len(arr)==4:
          atom = arr[0]
          x = float(arr[1])
          y = float(arr[2])
          z = float(arr[3])
          data.append((atom, x, y, z))
  return data

 
def read_central_and_blocks(atom_file):
  try:
    with open(atom_file, 'r') as fin:
      raw_input = fin.read().splitlines()
  except FileNotFoundError:
    print('file not found')
  center_number = raw_input[0].split(':')[1]
  blocks_dict = {}
  for block in raw_input[1:]:
    block_name = block.split(':')[0]
    block_atoms = block.split(':')[1].split(',')
    blocks_dict.setdefault(block_name, block_atoms)
  return center_number, blocks_dict


def atom_type_number(atom_xyz, index):
  return atom_xyz[index][0]+str(index+1)
  
  
def get_blocks_dists(atom_xyz, center_num, blocks):
  '''
  type center_num: str(int)
  type block: dict{str: list[str]}
  rtype: dict{str{str: float}}
  '''
  center_index = int(center_num) - 1
  center_atom = atom_type_number(atom_xyz, center_index)
  block_indexes = {k: [int(v)- 1 for v in vs] for k, vs in blocks.items()}
  block_pair_dists = {k : {center_atom + '_' + atom_type_number(atom_xyz, v) : dist(atom_xyz, center_index, v) for v in block_indexes[k]} for k in block_indexes.keys()}
  return block_pair_dists

def wrap_blocks_dists(comfiles, center_atom_num, blocks_dict):
  all_dists = []
  for com in comfiles:
    input = com
    atom_xyz = read_onefile(input)
    all_dists.append(get_blocks_dists(atom_xyz, center_atom_num, blocks_dict))
  return all_dists

def write_distance_file(comfiles, all_dists, outputname):
  with open(outputname, 'w') as fout:
    fout.write('\t'.join(['file','block','pair','distance'])+'\n')
    for i in list(range(len(comfiles))):
      com_name = comfiles[i]
      curr_dict = all_dists[i]
      for block, pair_dist in curr_dict.items():
        for pair, dist in pair_dist.items():
          fout.write('\t'.join([com_name, str(block), str(pair), str(dist)]) + '\n')
  
def filter_count(block_pair_dists):
  block_filter = {k : [dist <= FILTER for dist in block_pair_dists[k].values()] for k in block_pair_dists.keys()}
  block_sum = {k : sum(bool_list) for k, bool_list in block_filter.items()}
  total_count = sum(block_sum.values())
  
  if total_count == 0:
    block_sum.setdefault('TOTAL_COUNT',total_count)
    return block_sum
  block_frac = {k : v/total_count for k, v in block_sum.items()}
  block_frac.setdefault('TOTAL_COUNT', total_count)
  
  return block_frac


def write_fraction(df, all_dists, frac_fout):
  list_frac = list(map(filter_count, all_dists))
  frac_df =pd.DataFrame(list_frac)
  new_df = pd.concat([df.iloc[:,:-2], frac_df, df.iloc[:, -2:]], axis = 1)
  new_df.to_csv(frac_fout, sep='\t',index=False)

  
if __name__ == '__main__':
  if len(sys.argv) < 3:
    print('Usage: python filter_distance.py kmean_file atom_file')
    sys.exit()
  
  kmean_file = sys.argv[1]
  atom_file = sys.argv[2]
  center_atom_num, blocks_dict = read_central_and_blocks(atom_file)
  print('central atom is {}'.format(center_atom_num))
  for k, v in blocks_dict.items():
    print('{0} has {1}'.format(k, v))
  df = pd.read_csv(kmean_file, sep = '\t')
  comfiles = df['file']
  all_dists = wrap_blocks_dists(comfiles, center_atom_num, blocks_dict)
  dist_fout = 'block_pair_distances.tsv'
  str_filter = str(FILTER).replace('.','p')
  frac_fout = kmean_file[:-4] + '_block_frac_filter'+str_filter+'.tsv'
  write_distance_file(comfiles, all_dists, dist_fout)
  write_fraction(df, all_dists,frac_fout)
    

  
    
  
