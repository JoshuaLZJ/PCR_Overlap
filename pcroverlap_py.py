import primer3
import pandas as pd
from fasta_reader import read_fasta
import argparse

parser = argparse.ArgumentParser(description='Generate overlapping primers for a specified region')
parser.add_argument('-f', '--file', metavar='<path-to-fasta>', type=str, help='file path for fasta file')
parser.add_argument('-r', '--size_range', metavar='<int>', type=int, nargs=2, help='size range for PCR products')
parser.add_argument('-l', '--overlap', metavar='<int>', type=int, help='size of overlap between PCR products')
parser.add_argument('-s', '--start', metavar='<int>', type=int, help='start site for generating overlapping primers in fasta sequence')
parser.add_argument('-o', '--output_file', metavar='<path-to-csv>', type=str, help='file path for output')
parser.print_help()

args = parser.parse_args()

file = args.file
size_range = args.size_range
overlap = args.overlap
start = args.start
output_file = args.output_file

fasta = read_fasta(file)
item = fasta.read_item()

res_lst = []

anchor = start

global_arg = {
          'PRIMER_OPT_SIZE': 20,
          'PRIMER_PICK_INTERNAL_OLIGO': 1,
          'PRIMER_INTERNAL_MAX_SELF_END': 8,
          'PRIMER_MIN_SIZE': 18,
          'PRIMER_MAX_SIZE': 25,
          'PRIMER_OPT_TM': 60.0,
          'PRIMER_MIN_TM': 57.0,
          'PRIMER_MAX_TM': 63.0,
          'PRIMER_MIN_GC': 20.0,
          'PRIMER_MAX_GC': 80.0,
          'PRIMER_MAX_POLY_X': 100,
          'PRIMER_INTERNAL_MAX_POLY_X': 100,
          'PRIMER_SALT_MONOVALENT': 50.0,
          'PRIMER_DNA_CONC': 50.0,
          'PRIMER_MAX_NS_ACCEPTED': 0,
          'PRIMER_MAX_SELF_ANY': 12,
          'PRIMER_MAX_SELF_END': 8,
          'PRIMER_PAIR_MAX_COMPL_ANY': 12,
          'PRIMER_PAIR_MAX_COMPL_END': 8,
          'PRIMER_PRODUCT_SIZE_RANGE': [size_range],
      }
i=0
while (anchor + size_range[1]) < len(item.sequence):
  try:
    primer_search_regions = [anchor, overlap, anchor + size_range[1], overlap]
    print('Searching Primers in region: ', primer_search_regions[0], '-', primer_search_regions[0]+primer_search_regions[1], 'and', primer_search_regions[2], '-', primer_search_regions[2]+primer_search_regions[3])
    if i == 0:
      seq_arg = {
            'SEQUENCE_ID': 'testfa',
            'SEQUENCE_TEMPLATE': item.sequence,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': primer_search_regions,
        }
    else:
      seq_arg = {
            'SEQUENCE_ID': 'testfa',
            'SEQUENCE_TEMPLATE': item.sequence,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': primer_search_regions,
            'SEQUENCE_EXCLUDED_REGION': excluded_region
        }
    res = primer3.bindings.designPrimers(global_args=global_arg, seq_args=seq_arg)
    if res['PRIMER_LEFT_NUM_RETURNED'] > 0:
      res_lst.append(res)
    else:
      shift = 100
      while res['PRIMER_LEFT_NUM_RETURNED'] == 0:
        print('No primers discovered in the region')
        shift += 100
        primer_search_regions = [anchor - shift, overlap, anchor + size_range[1] - shift, overlap]
        print('Searching Primers in region: ', primer_search_regions[0], '-', primer_search_regions[0]+primer_search_regions[1], 'and', primer_search_regions[2], '-', primer_search_regions[2]+primer_search_regions[3])
        seq_arg = {
            'SEQUENCE_ID': 'testfa',
            'SEQUENCE_TEMPLATE': item.sequence,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': primer_search_regions,
            'SEQUENCE_EXCLUDED_REGION': excluded_region
        }
        res = primer3.bindings.designPrimers(global_args=global_arg, seq_args=seq_arg)
      res_lst.append(res)
    anchor = res['PRIMER_RIGHT_0'][0] - res['PRIMER_RIGHT_0'][1] - overlap
    excluded_region = [res['PRIMER_RIGHT_0'][0] -  res['PRIMER_RIGHT_0'][1], res['PRIMER_RIGHT_0'][1]]
    i += 1
  except:
    print('program failed')
    break

res_df = pd.DataFrame(res_lst)
res_df.to_csv(output_file)
print('All overlapping primers saved in: ', output_file)