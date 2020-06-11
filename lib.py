import numpy as np, pandas as pd
import plotly, os, sys
from collections import defaultdict

app_fold = os.path.dirname(os.path.realpath(__file__)) + '/'

## Parameters
efficiency_model_std = 1.5

###############################################
# Functions and variables for URL shortening
###############################################

chars = None
dna_to_code = dict()
code_to_dna = dict()

KMER_LEN = 9

def __init_chars():
  global chars
  chars = [chr(s) for s in range(48, 48 + 10)] + [chr(s) for s in range(65, 65 + 26)] + [chr(s) for s in range(97, 97 + 26)]
  chars += ['-', '_', '~', '.']
  chars.remove('_')
  return

def __init_mappers():
  output = chars
  # All 3-mers of 65-length safe html character alphabet
  for idx in range(3-1):
    output = __append_alphabet(output, chars)
  triplets = output

  # All 9-mers DNA
  output = list('ACGT')
  for idx in range(KMER_LEN-1):
    output = __append_alphabet(output, list('ACGT'))
  kmers = output

  global dna_to_code
  global code_to_dna
  for kmer, triplet in zip(kmers, triplets):
    dna_to_code[kmer] = triplet
    code_to_dna[triplet] = kmer
  return

def __append_alphabet(output, alphabet):
  new_output = []
  for o in output:
    for a in alphabet:
      new_output.append(o + a)
  return new_output

def parse_coded_seq_leftover(encoded_dna, leftover_dna):
  # Process encoded DNA
  if len(leftover_dna) != 1 and len(encoded_dna) % 3 != 0:
    return '-'
  if encoded_dna == '-':
    return leftover_dna

  seq = ''
  for jdx in range(0, len(encoded_dna), 3):
    w = encoded_dna[jdx : jdx + 3]
    seq += code_to_dna[w]
  if leftover_dna != '-':
    seq += leftover_dna
  return seq

def encode_dna(seq):
  if seq is None or len(seq) == 0:
    return '-', '-'

  if len(seq) < KMER_LEN:
    return '-', seq

  encodeddna = ''
  for idx in range(0, len(seq), KMER_LEN):
    chomp = seq[idx : idx + KMER_LEN]
    if len(chomp) == KMER_LEN:
      encodeddna += dna_to_code[chomp]
    else:
      break
  if len(seq[idx:]) != KMER_LEN:
    leftoverdna = seq[idx:]
  else:
    leftoverdna = '-'
  return encodeddna, leftoverdna

###############################################
# Single
###############################################

def parse_valid_url_path_single(url_path):
  ## Expected format:
  if url_path[:len('/s_')] != '/s_': return False, None

  url_path = url_path.replace('/s_', '')
  if len(url_path) == 0 or '_' not in url_path: return False, None

  parts = url_path.split('_')
  cats = [
    'base_editor', 
    'celltype', 
    'encodeddna', 
    'leftoverdna', 
    'psidx', 
    'aa_frame_txt'
  ]

  if len(parts) != len(cats): return False, None
  results_d = {cat: parts[idx] for idx, cat in enumerate(cats)}

  seq = parse_coded_seq_leftover(
    results_d['encodeddna'], 
    results_d['leftoverdna'],
  )
  results_d['seq'] = seq

  psidx = int(results_d['psidx'])
  results_d['protospacer'] = f'{seq[psidx : psidx + 20]}, {psidx}'

  aaft = results_d['aa_frame_txt']
  frame = int(aaft[0])
  mapper = {
    'p': '+',
    'm': '-',
  }
  strand = mapper[aaft[1]]
  results_d['aa_frame_txt_parsed'] = f'{frame},{strand}'

  return True, results_d

def encode_dna_to_url_path_single(seq, ps, base_editor, celltype, aa_frame_txt):
  seq = seq.upper()
  encodeddna, leftoverdna = encode_dna(seq)
  [protospacer, ps_idx] = ps.split()

  if aa_frame_txt != 'None':
    aa_frame = int(aa_frame_txt[0])
    aa_strand = aa_frame_txt[-1]
    strand_to_text = {
      '+': 'p',
      '-': 'm',
    }
    aa_frame_txt = f'{aa_frame}{strand_to_text[aa_strand]}'

  return f'/s_{base_editor}_{celltype}_{encodeddna}_{leftoverdna}_{ps_idx}_{aa_frame_txt}'


###############################################
# Batch
###############################################

# def parse_valid_url_path_batch(url_path):
#   ## Expected format:
#   # [encodedDNA]_[leftoverDNA]_[pam in plaintext] + more
#   dd = dict()
#   if url_path[:len('/batch_')] != '/batch_':
#     return False, dd

#   url_path = url_path.replace('/batch_', '')
#   if len(url_path) == 0 or '_' not in url_path:
#     return False, dd

#   parts = url_path.split('_')
#   cats = ['celltype', 'coded', 'leftover', 'pam', 'adv_flag', 'coded_spec', 'leftover_spec', 'adv_poi', 'adv_delstart', 'adv_delend', 'chosen_columns', 'sort_by', 'sort_dir', 'row_select']
#   if len(parts) != len(cats):
#     return False, dd
#   for idx, cat in enumerate(cats):
#     dd[cat] = parts[idx]

#   dd['seq'] = parse_coded_seq_leftover(dd, 'coded', 'leftover')
#   dd['adv_seq_spec'] = parse_coded_seq_leftover(dd, 'coded_spec', 'leftover_spec')

#   # Reword some values
#   if dd['adv_flag'] == '1':
#     dd['adv_flag'] = True
#   elif dd['adv_flag'] == '0':
#     dd['adv_flag'] = False

#   if dd['sort_dir'] == '1':
#     dd['sort_dir'] = 'Ascending'
#   else:
#     dd['sort_dir'] = 'Descending'

#   return True, dd

# def encode_dna_to_url_path_batch(seq, pam, celltype, adv_flag, adv_seq_spec, adv_poi, adv_delstart, adv_delend, chosen_columns, column_options, sort_by, sort_dir, selected_row):
#   seq, pam = seq.upper(), pam.upper()
#   edna, ldna = encode_dna(seq)
#   edna2, ldna2 = encode_dna(adv_seq_spec)

#   if adv_flag == True:
#     adv_flag_val = '1'
#   else:
#     adv_flag_val = '0'

#   adv_poi = transform_empty_value_to_dash(adv_poi)
#   adv_delstart = transform_empty_value_to_dash(adv_delstart)
#   adv_delend = transform_empty_value_to_dash(adv_delend)
#   sort_by = transform_empty_value_to_dash(sort_by)

#   binary_flags_chosen_cols = ''
#   for co in sorted([s['value'] for s in column_options]):
#     if co in chosen_columns:
#       binary_flags_chosen_cols += '1'
#     else:
#       binary_flags_chosen_cols += '0'

#   if sort_by != '-':
#     sort_by = sorted(chosen_columns).index(sort_by)

#   if sort_dir == 'Ascending':
#     sort_dir_val = '1'
#   else:
#     sort_dir_val = '0'

#   if selected_row == []:
#     selected_row_val = '-'
#   else:
#     selected_row_val = selected_row[0]

#   items = [
#     celltype,
#     edna, 
#     ldna, 
#     pam, 
#     adv_flag_val, 
#     edna2, 
#     ldna2, 
#     adv_poi, 
#     adv_delstart, 
#     adv_delend, 
#     binary_flags_chosen_cols, 
#     sort_by, 
#     sort_dir_val, 
#     selected_row_val
#   ]
#   return '/batch_%s' % ('_'.join([str(s) for s in items]))

def transform_empty_value_to_dash(val):
  if val is None or len(val) == 0 or val == 'None':
    return '-'
  else:
    return val

__init_chars()
__init_mappers()


###############################################
# Gene
###############################################

def parse_valid_url_path_gene(url_path):
  dd = dict()
  if url_path[:len('/gene_')] != '/gene_':
    return False, dd

  url_path = url_path.replace('/gene_', '')
  if len(url_path) == 0 or '_' not in url_path:
    return False, dd

  parts = url_path.split('_')
  cats = ['genome_build', 'gene', 'celltype', 'chosen_columns', 'sort_by', 'sort_dir', 'row_select']
  if len(parts) != len(cats):
    return False, dd
  for idx, cat in enumerate(cats):
    dd[cat] = parts[idx]

  if dd['sort_dir'] == '1':
    dd['sort_dir'] = 'Ascending'
  else:
    dd['sort_dir'] = 'Descending'

  return True, dd


def encode_url_path_gene(genome_build, gene, celltype, chosen_columns, column_options, sort_by, sort_dir, selected_row):
  binary_flags_chosen_cols = ''
  for co in sorted([s['value'] for s in column_options]):
    if co in chosen_columns:
      binary_flags_chosen_cols += '1'
    else:
      binary_flags_chosen_cols += '0'

  if sort_by is not None:
    sort_by = sorted(chosen_columns).index(sort_by)
  else: 
    sort_by = '-'

  if sort_dir == 'Ascending':
    sort_dir_val = '1'
  else:
    sort_dir_val = '0'

  if selected_row == []:
    selected_row_val = '-'
  else:
    selected_row_val = selected_row[0]

  items = [
    genome_build,
    gene,
    celltype,
    binary_flags_chosen_cols, 
    sort_by, 
    sort_dir_val, 
    selected_row_val,
  ]
  return '/gene_%s' % ('_'.join([str(s) for s in items]))

###############################################
# Data management
###############################################
bystander_models_design = pd.read_csv(app_fold + f'be_predict_bystander/models.csv', index_col = 0)
editor_celltype_dropdown_options = []
batch_dropdown_options = []
type_to_bes = defaultdict(list)
'''
  options = [
    {'label': '1-bp insertions', 'value': '1-bp insertions'},
'''
for idx, row in bystander_models_design.iterrows():
  if row['Hidden']:
    continue
  editor = row['Public base editor']
  celltype = row['Celltype']
  be_type = row['Public base editor type']
  editor_celltype_dropdown_options.append({
    'label': f'{editor}, {celltype}',
    'value': f'{editor}, {celltype}',
  })

  batch_option = {
    'label': f'{be_type}, {celltype}',
    'value': f'{be_type}, {celltype}',
  }
  if batch_option not in batch_dropdown_options:
    batch_dropdown_options.append(batch_option)

  type_to_bes[f'{be_type}, {celltype}'].append(editor)

###############################################
# Compbio operations
###############################################

def revcomp(seq):
  rc_mapper = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
  rc_seq = []
  for c in seq:
    if c in rc_mapper:
      rc_seq.append(rc_mapper[c])
    else:
      rc_seq.append(c)
  return ''.join(rc_seq[::-1])

def pam_shift(text1, text2, text_pam, direction):
  seq = text1 + text2
  cutsite = len(text1)

  if direction == 'right':
    cutsites = range(cutsite + 1, len(seq))
  elif direction == 'left':
    cutsites = range(cutsite - 1, 0, -1)

  for ct in cutsites:
    candidate_pam = seq[ct + 3 : ct + 6]
    if match(text_pam, candidate_pam):
      return seq[:ct], seq[ct:]
  return None

mapper = {
  'A': list('A'),
  'C': list('C'),
  'G': list('G'),
  'T': list('T'),
  'Y': list('CT'),
  'R': list('AG'),
  'W': list('AT'),
  'S': list('GC'),
  'K': list('TG'),
  'M': list('AC'),
  'D': list('AGT'),
  'V': list('ACG'),
  'H': list('ACT'),
  'B': list('CGT'),
  'N': list('ACGT'),
}
def match(template, dna):
  if len(dna) != len(template):
    return False
  for char, t in zip(dna, template):
    if char not in mapper[t]:
      return False
  return True

def estimate_pam_freq(pam):
  factor = 1
  for char in pam:
    factor *= ( len(mapper[char]) / 4)
  return factor

# codon_table = {
#   'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
#   'Leu': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
#   'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
#   'Lys': ['AAA', 'AAG'],
#   'Asn': ['AAT', 'AAC'], 
#   'Met': ['ATG'],
#   'Asp': ['GAT', 'GAC'], 
#   'Phe': ['TTT', 'TTC'],
#   'Cys': ['TGT', 'TGC'],
#   'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
#   'Gln': ['CAA', 'CAG'],
#   'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
#   'Glu': ['GAA', 'GAG'],
#   'Thr': ['ACT', 'ACC', 'ACA', 'ACG'],
#   'Gly': ['GGT', 'GGC', 'GGA', 'GGG'],
#   'Trp': ['TGG'],
#   'His': ['CAT', 'CAC'], 
#   'Tyr': ['TAT', 'TAC'],
#   'Ile': ['ATT', 'ATC', 'ATA'],
#   'Val': ['GTT', 'GTC', 'GTA', 'GTG'],
#   'Stop': ['TAA', 'TGA', 'TAG'],
# }
codon_table = {
  'A': ['GCT', 'GCC', 'GCA', 'GCG'],
  'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
  'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
  'K': ['AAA', 'AAG'],
  'N': ['AAT', 'AAC'], 
  'M': ['ATG'],
  'D': ['GAT', 'GAC'], 
  'F': ['TTT', 'TTC'],
  'C': ['TGT', 'TGC'],
  'P': ['CCT', 'CCC', 'CCA', 'CCG'],
  'Q': ['CAA', 'CAG'],
  'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
  'E': ['GAA', 'GAG'],
  'T': ['ACT', 'ACC', 'ACA', 'ACG'],
  'G': ['GGT', 'GGC', 'GGA', 'GGG'],
  'W': ['TGG'],
  'H': ['CAT', 'CAC'], 
  'Y': ['TAT', 'TAC'],
  'I': ['ATT', 'ATC', 'ATA'],
  'V': ['GTT', 'GTC', 'GTA', 'GTG'],
  '*': ['TAA', 'TGA', 'TAG'],
}
def dna_to_aa(dna, frame, strand):
  '''
    frame in [0, 1, 2]
    strand in ['+', '-']
  '''
  dna = dna.upper()
  aas = ''
  if strand == '-': 
    dna = revcomp(dna)

  start_gap = frame
  dna = dna[frame:]
  end_gap = 0
  for idx in range(0, len(dna), 3):
    codon = dna[idx : idx + 3]
    if len(codon) != 3: 
      end_gap = len(codon)
      break
    for aa in codon_table:
      if codon in codon_table[aa]:
        aas += aa
        break

  if strand == '-':
    aas = aas[::-1]
  return aas

def get_aa_display_start_idx(dna, frame, strand):
  '''
    frame in [0, 1, 2]
    strand in ['+', '-']
  '''
  dna = dna.upper()
  aas = ''
  if strand == '-': 
    dna = revcomp(dna)

  start_gap = frame
  dna = dna[frame:]
  end_gap = 0
  for idx in range(0, len(dna), 3):
    codon = dna[idx : idx + 3]
    if len(codon) != 3: 
      end_gap = len(codon)
      break

  if strand == '-':
    start_gap, end_gap = end_gap, start_gap
  return start_gap, end_gap

def find_protospacers(seq):
  left_margin = 20
  right_margin = 10
  ps_len = 20
  protospacers = []
  poss = []
  unique_flags = []
  for idx in range(left_margin, len(seq) - ps_len - right_margin):
    ps = seq[idx : idx + ps_len]
    unique_flags.append(bool(ps not in protospacers))
    protospacers.append(ps)
    poss.append(idx)
  return protospacers, poss, unique_flags

###############################################
# Colors
###############################################
rgb = {
  'red': 'rgb(236, 67, 57)',
  'orange': 'rgb(244, 123, 22)',
  'light green': 'rgb(174, 214, 119)',
  'green': 'rgb(124, 184, 47)',
  'dark green': 'rgb(78, 143, 19)',
  'blue': 'rgb(0, 174, 179)',
  'dark blue': 'rgb(0, 136, 145)',
  'lilac': 'rgb(191, 171, 230)',
  'pink': 'rgb(243, 113, 175)',
  'gray': 'rgb(134, 137, 140)',
}

editor_cmap = {
  'ABE': rgb['red'],
  'ABE8': rgb['red'],
  'BE4': 'rgb(0, 160, 220)',
  'BE4-CP1028': rgb['blue'],
  'eA3A': 'rgb(140, 104, 203)',
  'AID': rgb['green'],
  'CDA': 'rgb(239, 185, 23)',
  'ABE-CP1041': rgb['orange'],
  'evoAPOBEC': 'rgb(237, 71, 149)',
}

aa_cmap = {
  'A': 'rgba(174, 214, 119, 0.5)',
  'G': 'rgba(199, 229, 154, 0.5)',
  'C': 'rgba(149, 199, 83, 0.5)',
  'N': 'rgba(124, 184, 47, 0.5)',
  'D': 'rgba(96, 170, 20, 0.5)',
  'Q': 'rgba(78, 143, 19, 0.5)',
  'E': 'rgba(59, 117, 71, 0.5)',
  'S': 'rgba(236, 67, 57, 0.5)',
  'T': 'rgba(241, 109, 100, 0.5)',
  'R': 'rgba(244, 123, 22, 0.5)',
  'K': 'rgba(245, 150, 64, 0.5)',
  'M': 'rgba(134, 137, 140, 0.5)',
  'L': 'rgba(0, 174, 179, 0.5)',
  'I': 'rgba(53, 190, 193, 0.5)',
  'V': 'rgba(0, 160, 220, 0.5)',
  'H': 'rgba(0, 140, 201, 0.5)',
  'F': 'rgba(216, 204, 244, 0.5)',
  'Y': 'rgba(191, 171, 230, 0.5)',
  'W': 'rgba(165, 137, 217, 0.5)',
  'P': 'rgba(243, 113, 175, 0.5)',
  '*': 'rgba(115, 118, 121, 0.5)',
  ' ': 'white',
}

# Used for making color scale for [0, 1] values
dna_color_minmax = {
  'A': [
    'rgba(255, 224, 218, 1)',
    'rgba(221, 46, 31, 1)',
  ],
  'C': [
    'rgba(255, 242, 182, 1)',
    'rgba(230, 167, 0, 1)',
  ],
  'G': [
    'rgba(224, 244, 190, 1)',
    'rgba(96, 170, 20, 1)',
  ],
  'T': [
    'rgba(207, 237, 251, 1)',
    'rgba(0, 140, 201, 1)',
  ],
}

# Font color for match vs. edited nts / aas
font_cmap = {
  'match': 'rgba(208, 211, 214, 1)',
  'edited': 'black',
}

num_colors_in_ref = 666
num_colors_in_ref_resolution = 1 / num_colors_in_ref
dna_color_scales = {
  nt: plotly.colors.n_colors(
    dna_color_minmax[nt][0],
    dna_color_minmax[nt][1],
    num_colors_in_ref, 
    colortype = 'rgb'
  ) for nt in dna_color_minmax
}

gray_color_scale = plotly.colors.n_colors(
  'rgba(230, 233, 236, 1)',
  'rgba(115, 118, 121, 1)',
  num_colors_in_ref, 
  colortype = 'rgb'
)

def get_color(scale, val, white_threshold = 0):
  # val in [0, 1]. scale = list of colors. assumed to be more white at 0
  if val < white_threshold: return 'white'
  c_idx = int(val / num_colors_in_ref_resolution)    
  return scale[c_idx]

###############################################
# 
###############################################

###############################################
# Alignment text presentation
###############################################
'''
  Keep to 4 chars max
'''
editor_to_header = {
  'BE4': ['BE4'],
  'CDA': ['CDA'],
  'AID': ['AID'],
  'evoAPOBEC': ['evoA'],
  'ABE': ['ABE'],
  'eA3A': ['eA3A'],
  'eA3A-T31A': ['eA3A', 'T31A'],
  'eA3A-T31AT44A': ['eA3A', 'T31A', 'T44A'],
  'eA3A-BE5': ['eA3A', 'BE5'],
  'BE4-H47ES48A': ['BE4', 'H47E', 'S48A'],
  'EA-BE4': ['EA-', 'BE4'],
  'ABE8': ['ABE', '8e'],
  'ABE-CP1041': ['ABE-', 'CP', '1041'],
  'BE4-CP1028': ['BE4-', 'CP', '1028'],
}


###############################################
# Batch mode: xaxis ticks 
###############################################
def get_batch_statcol_xrange(stats, stat_nm):
  if '(%)' in stat_nm:
    buff = 3
  elif stat_nm in ['Exp. indel len', 'Exon number']:
    buff = 1
  elif stat_nm == 'MH strength':
    buff = 0.1
  elif stat_nm == 'Precision':
    buff = 0.05 
  elif stat_nm in ['Cutsite', 'Dist. to 5\' end', 'Dist. to 3\' end']:
    buff = 10
  elif stat_nm in ['Repairs to spec.', 'Deletes spec.']:
    buff = 5
  elif stat_nm == 'Dist. to POI':
    buff = 5
  else: # default
    buff = 0
  return [min(stats) - buff, max(stats) + buff]

# def get_batch_statcol_xticks(stats):
  # pass
  # return

def get_batch_select_line(x0 = 0, x1 = 0, y0 = 0, y1 = 0, xref = '', yref = ''):
  return dict(
    type = 'line',
    xref = xref,
    yref = yref,
    x0 = x0,
    x1 = x1,
    y0 = y0,
    y1 = y1,
    opacity = 0.8,
    line = dict(
      color = 'rgb(33, 33, 33)',
      width = 1,
      dash = 'dot',
    )
  )

def rename_batch_columns(stats):
  name_changes = {
    'Frameshift frequency': 'Frameshift (%)',
    'Frame +0 frequency': 'Frame +0 (%)',
    'Frame +1 frequency': 'Frame +1 (%)',
    'Frame +2 frequency': 'Frame +2 (%)',
    'Highest outcome frequency': 'M.F. gt (%)',
    'Highest del frequency': 'M.F. del (%)',
    'Highest ins frequency': 'M.F. ins (%)',
    'Expected indel length': 'Exp. indel len',
    'Distance to 5\' exon boundary': 'Dist. to 5\' end',
    'Distance to 3\' exon boundary': 'Dist. to 3\' end',
  }
  for col in stats:
    if col in name_changes:
      stats[name_changes[col]] = stats[col]
      stats.drop([col], axis = 1, inplace = True)
  return stats

def order_chosen_columns(cols):
  preferred_order = [
    'Exon number',
    'Dist. to 5\' end',
    'Dist. to 3\' end',
    'Cutsite',
    'Dist. to POI',
    'Repairs to spec.',
    'Deletes spec.',
    'Precision',
    'Frameshift (%)',
    'Frame +0 (%)',
    'Frame +1 (%)',
    'Frame +2 (%)',
    'MH strength',
    'M.F. gt (%)',
    'M.F. del (%)',
    'M.F. ins (%)',
    'Exp. indel len',
  ]
  reordered = []
  for pref in preferred_order:
    if pref in cols:
      reordered.append(pref)
  return reordered

def get_x_domains(num_cols):
  # Ensure uniform and consistent horizontal spacing with variable number of columns
  margin_pct = 0.12

  domains = []
  for leftside in np.arange(0, 1, 1/num_cols):
    size = 1 / num_cols
    margin_size = size * margin_pct 
    rightside = leftside + size
    domains.append([leftside + margin_size, rightside - margin_size])
  return domains

def get_fixedwidth_ID(ids):
  largest_len = len(str(max(ids)))
  fw_ids = []
  for item in ids:
    num_spaces = largest_len - len(str(item))
    fw_id = '%s#%s' % (' ' * num_spaces, item)
    fw_ids.append(fw_id)
  return fw_ids

def get_fixedwidth_items(items):
  largest_len = len(str(max(items)))
  fw_items = []
  for item in items:
    num_spaces = largest_len - len(str(item))
    fw_item = '%s%s' % (' ' * num_spaces, item)
    fw_items.append(fw_item)
  return fw_items