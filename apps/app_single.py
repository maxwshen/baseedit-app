import pickle, copy, os, datetime, subprocess, json, sys
from collections import defaultdict
import random
import numpy as np
import pandas as pd
from scipy.stats import entropy
import time
from io import StringIO
from textwrap import dedent

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import flask
import plotly
from flask_caching import Cache

import lib, header

from app_holder import app

# Import models
app_fold = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/'
sys.path.append(app_fold)
from be_predict_bystander import predict as bystander_model
from be_predict_efficiency import predict as efficiency_model

bystander_model.init_all_models()
efficiency_model.init_all_models()

try:
  os.mkdir('user-csvs/')
except FileExistsError:
  pass
else:
  subprocess.check_output('rm -rf user-csvs/*', shell = True)

# Set up flask caching
CACHE_CONFIG = {
  'CACHE_TYPE': 'redis',
  'CACHE_REDIS_URL': os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()
cache.init_app(app.server, config = CACHE_CONFIG)
cache_timeout = 300

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

# Random default, which is cached on filesystem
default_text = ''.join([random.choice(list('ACGT')) for s in range(50)])
if os.path.isfile('single_default.pkl'):
  subprocess.check_output('rm -rf single_default.pkl', shell = True)


## Parameters

###################################################################
###################################################################
##
# App layout
##
layout = html.Div([

  ##
  # Hidden divs for light data storage
  ##
  html.Div(
    [
      html.Div(
        id = 'S_hidden-pred-signal_bystander',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden-pred-signal_efficiency',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden-chosen-base_editor',
        children = 'BE4',
      ),
      html.Div(
        id = 'S_hidden-chosen-celltype',
        children = 'mES',
      ),

      dcc.Location(
        id = 'S_url',
        refresh = False,
      ),
    ],
    style = dict(
      display = 'none',
    ),
  ),

  ##
  # Header
  ##
  html.Div(
    [
      ###################################################
      # Upper header
      ###################################################
      header.get_navigation_header('single'),

      ###################################################
      # Sequence boxes
      ###################################################
      html.Div([
        html.Div(
          [
            dcc.Input(
              id = 'S_textbox', 
              size = '28',
              value = default_text,
              type = 'text',
              autoFocus = True,
              style = dict(
                fontFamily = 'monospace',
                fontSize = 16,
                # float = 'left',
                transform = 'translateX(50%)',
              ),
            )
          ],
          className = 'dna_textbox',
        ),
        ], style = dict(
          verticalAlign = 'center',
          whiteSpace = 'nowrap',
          overflowX = 'auto',
        ),
      ),

      # Empty div for bottom margin in header
      html.Div(
        [], 
        style = dict(
          marginBottom = '10px',
        ),
      )

    ],
    style = dict(
      position = 'fixed',
      top = 0,
      backgroundColor = 'white',
      borderBottom = '3px solid #777777',
      zIndex = 1e6,
      width = '1010px',
      left = '50%',
      transform = 'translate(-50%, 0px)',
    ),
  ),

  ##
  # Body / plots
  ##
  html.Div(
    [
      # First (not bottom animated)
      ###################################################
      # Module: Single base editor
      ###################################################
      html.Div([
        # header
        html.Div([
          html.Div([
            html.Strong('Ex: Single base editor')
            ],
            className = 'module_header_text'),
          ],
          className = 'module_header'
        ),

        html.Div(
          id = 'S_text-test_efficiency',
        ),

        html.Div(
          [
            # Text table
            dcc.Graph(
              id = 'S_summary-alignment-table',
              config = dict(
                modeBarButtonsToRemove = modebarbuttons_2d,
                displaylogo = False,
                displayModeBar = False,
                staticPlot = True,
              ),
              style = dict(
                height = 290,
                width = 629,
              ),
              className = 'twelve columns',
            ),

          ],
          className = 'row',
        ),

      ], className = 'module_style',
      ),

      # Animate bottom
      html.Div([

        ###################################################
        # Module: Single base editor
        ###################################################
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('Ex: Single base editor')
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),
        ], className = 'module_style',
        ),

        ],
        id = 'S_plots_body',
        style = dict(
          display = 'none',
        ),
        className = 'animate-bottom',
      ),

    ],
    # body style
    # id = 'S_plots_body',
    style = dict(
      # display = 'none',
      transform = 'translateY(%spx)' % (200),
    ),
  ),
  ##

  ],  # body div
  style = dict(
    width = '970px',
    margin = '0 auto',
  )
)

#######################################################################
#########################      CALLBACKS      #########################
#######################################################################

##
# Celltype choice callbacks
## 

##
# Prediction caching
##
@cache.memoize(timeout = cache_timeout)
def bystander_predict_cache(seq, base_editor, celltype):
  pred_df, stats = bystander_model.predict_given(
    seq,
    base_editor = base_editor,
    celltype = celltype,
  )
  pred_df = bystander_model.add_genotype_column(pred_df, stats)
  filtered_cols = ['Predicted frequency', 'Genotype']
  nt_cols = [col for col in pred_df if col not in filtered_cols]
  return pred_df, stats, nt_cols

@cache.memoize(timeout = cache_timeout)
def efficiency_predict_cache(seq, base_editor, celltype):
  pred_d = efficiency_model.predict_given(
    seq,
    base_editor = base_editor,
    celltype = celltype,
  )
  return pred_d

##
# Prediction callbacks
##
@app.callback(
  Output('S_hidden-pred-signal_bystander', 'children'),
  [Input('S_textbox', 'value'),
   Input('S_hidden-chosen-base_editor', 'children'),
   Input('S_hidden-chosen-celltype', 'children')])
def bystander_predict(seq, base_editor, celltype):
  seq = seq.upper()
  bystander_predict_cache(seq, base_editor, celltype)
  return '%s,%s,%s' % (seq, base_editor, celltype)

@app.callback(
  Output('S_hidden-pred-signal_efficiency', 'children'),
  [Input('S_textbox', 'value'),
   Input('S_hidden-chosen-base_editor', 'children'),
   Input('S_hidden-chosen-celltype', 'children')])
def efficiency_predict(seq, base_editor, celltype):
  seq = seq.upper()
  efficiency_predict_cache(seq, base_editor, celltype)
  return '%s,%s,%s' % (seq, base_editor, celltype)

##
# Summary of predictions callbacks
##
@app.callback(
  Output('S_summary-alignment-table', 'figure'),
  [Input('S_hidden-pred-signal_bystander', 'children'),
  ])
def update_summary_alignment_text(signal):
  seq, base_editor, celltype = signal.split(',')
  pred_df, stats, nt_cols = bystander_predict_cache(seq, base_editor, celltype)

  ## Set up data
  top10 = pred_df.iloc[:10]
  fqs = top10['Predicted frequency']
  fq_strings = [''] + [f'{100*s:.1f}%' for s in fqs]

  p0idx = 19
  target_seq = stats['50-nt target sequence']
  target_aas = list(lib.dna_to_aa(target_seq, 0, '+')) + [' ']

  pred_df['Amino acid sequence'] = [lib.dna_to_aa(s, 0, '+') for s in pred_df['Genotype']]

  poswise_total = {col: sum(pred_df.loc[pred_df[col] != col[0], 'Predicted frequency']) for col in nt_cols}

  ## Set up colors
  fill_cmap = {
    'match': 'rgba(255, 255, 255, 1)',
    'A': 'rgba(236, 67, 57, 1)',
    'C': 'rgba(239, 185, 32, 1)',
    'G': 'rgba(124, 184, 47, 1)',
    'T': 'rgba(0, 160, 220, 1)',
  }

  font_cmap = {
    'match': 'rgba(208, 211, 214, 1)',
    'edited': 'black',
  }

  color_minmax = {
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

  num_colors_in_ref = 666
  num_colors_in_ref_resolution = 1 / num_colors_in_ref
  color_scales = {
    nt: plotly.colors.n_colors(
      color_minmax[nt][0],
      color_minmax[nt][1],
      num_colors_in_ref, 
      colortype = 'rgb'
    ) for nt in color_minmax
  }

  def get_color(scale, val, white_threshold = 0):
    # val in [0, 1]
    if val < white_threshold: return 'white'
    c_idx = int(val / num_colors_in_ref_resolution)    
    return scale[c_idx]

  ## Form table with colors
  fillcolors = []
  fontcolors = []
  poswise_cols = []
  for gt_idx, ref_nt in enumerate(target_seq):
    aa_idx = gt_idx // 3
    pos = gt_idx - p0idx
    cand_col = f'{nt}{pos}'
    pos_col = []
    col_fill_colors = []
    col_font_colors = []

    # Testing amino acid seq
    # Text
    ref_aa = target_aas[aa_idx]
    if gt_idx % 3 == 1:
      pos_col.append(ref_aa)
    else:
      pos_col.append('')
    # Color
    # col_fill_colors.append('rgba(0, 160, 220, 0)')
    col_fill_colors.append(lib.aa_cmap[ref_aa])
    col_font_colors.append('black')

    # row for target_seq 
    # Text
    pos_col.append(ref_nt)
    # Color
    if cand_col not in nt_cols:
      col_fill_colors.append('white')
    else:
      # col_fill_colors = [fill_cmap[ref_nt]]
      tot_edit_frac = poswise_total[cand_col]
      color_scale = color_scales[ref_nt]
      col_fill_colors.append(get_color(color_scale, tot_edit_frac, white_threshold = 0.0015))
    col_font_colors.append('black')

    # rows for edited genotypes + aas
    for jdx, row in top10.iterrows():
      pred_fq = row['Predicted frequency']
      gt_seq = row['Genotype']
      aa_seq = row['Amino acid sequence']

      obs_nt = gt_seq[gt_idx]
      obs_aa = aa_seq[aa_idx]

      # Amino acid row
      # Text
      if gt_idx % 3 == 1:
        pos_col.append(obs_aa)
      else: 
        pos_col.append('')
      # Color
      if obs_aa == ref_aa:
        col_fill_colors.append(fill_cmap['match'])
        col_font_colors.append(font_cmap['match'])
      else:
        col_fill_colors.append(lib.aa_cmap[obs_aa])
        col_font_colors.append(font_cmap['edited'])

      # Genotype row
      # Text
      # Color
      pos_col.append(obs_nt)
      if obs_nt == ref_nt:
        col_fill_colors.append(fill_cmap['match'])
        col_font_colors.append(font_cmap['match'])
      else:
        # col_fill_colors.append(fill_cmap[obs_nt])
        color_scale = color_scales[obs_nt]
        col_fill_colors.append(get_color(color_scale, pred_fq))
        col_font_colors.append(font_cmap['edited'])


    # Finished iterating over one column
    poswise_cols.append(pos_col)
    fillcolors.append(col_fill_colors)
    fontcolors.append(col_font_colors)

  # alignment_col_width = 420
  # pos_col_width = alignment_col_width // len(poswise_cols)
  pos_col_width = 1.2

  return dict(
    data = [go.Table(
      columnwidth = [pos_col_width] * len(poswise_cols) + [10],
      header = dict(
        line = dict(width = 0),
        fill = dict(color = 'white'),
        height = 0,
      ),
      cells = dict(
        values = poswise_cols + [fq_strings],
        align = ['center'] * len(poswise_cols) + ['right'], 
        fill = dict(
          color = fillcolors + ['rgba(255, 255, 255, 1)'] * len(fq_strings),
        ),
        line = dict(width = 0),
        font = dict(
          family = 'monospace',
          color = fontcolors + ['black'] * len(fq_strings),
        ),
        height = 20,
      ),
    )],
    layout = go.Layout(
      font = dict(
        family = 'monospace',
      ),
      margin = dict(
        l = 10,
        r = 0,
        t = 5,
        b = 5,
      ),
    ),
  )

##
# General stats callbacks
##

## General stats text
@app.callback(
  Output('S_text-test_efficiency', 'children'),
  [Input('S_hidden-pred-signal_efficiency', 'children')])
def text_genstats_precision(signal):
  seq, base_editor, celltype = signal.split(',')
  pred_d = efficiency_predict_cache(seq, base_editor, celltype)
  print(pred_d)
  logit_score = pred_d['Predicted logit score']

  return [
    html.Span(f'Logit score: {logit_score}'),
  ]

##
# Genotype table v2 callbacks
##

##
# Download callbacks
##

##
# Flask serving
##

##
# Page link callback
##
