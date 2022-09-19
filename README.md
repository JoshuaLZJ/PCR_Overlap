# PCR_Overlap
A Python script for designing overlapping primers to amplify a big genomic region.

Generate overlapping primers for a specified region

optional arguments:

-h, --help show this help message and exit
-f <path-to-fasta>, --file <path-to-fasta>, file path for fasta file
-r <int> <int>, --size_range <int> <int>, size range for PCR products
-l <int>, --overlap <int>, size of overlap between PCR products
-s <int>, --start <int>, start site for generating overlapping primers in fasta sequence
-o <path-to-csv>, --output_file <path-to-csv>, file path for output

usage: pcroverlap.py [-h] [-f <path-to-fasta>] [-r <int> <int>] [-l <int>] [-s <int>] [-o <path-to-csv>]
