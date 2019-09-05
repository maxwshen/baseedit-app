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
  bystander_model.init_model(
    base_editor = base_editor,
    celltype = celltype,
  )
  pred_df, stats = bystander_model.predict(seq)
  return pred_df, stats

@cache.memoize(timeout = cache_timeout)
def efficiency_predict_cache(seq, base_editor, celltype):
  efficiency_model.init_model(
    base_editor = base_editor,
    celltype = celltype,
  )
  pred_d = efficiency_model.predict(seq)
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
  pred_df, stats = bystander_predict_cache(seq, base_editor, celltype)

  top10 = pred_df.iloc[:10]
  fqs = top10['Predicted frequency']
  fq_strings = [''] + [f'{100*s:.1f}' for s in fqs]

  gt_cols = [col for col in pred_df.columns if col != 'Predicted frequency']
  p0idx = 19
  target_seq = stats['50-nt target sequence']
  alignments = [target_seq]
  print(top10.columns)
  for jdx, row in top10.iterrows():
    gt = ''
    for idx, nt in enumerate(target_seq):
      pos = idx - p0idx
      cand_col = f'{nt}{pos}'
      if cand_col not in gt_cols:
        # bullet: •. Use with smaller font size
        gt += '.'
      else:
        gt += row[cand_col]
    alignments.append(gt)

  return dict(
    data = [go.Table(
      columnwidth = [310, 110, 60],
      header = dict(
        values = ['Alignment', '%'],
        align = ['center', 'right'], 
        line = dict(color = 'white'),
        fill = dict(color = 'white'),
      ),
      cells = dict(
        values = [alignments, fq_strings],
        align = ['center', 'right'], 
        line = dict(color = 'white'),
        fill = dict(color = 'white'),
        font = dict(family = 'monospace'),
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
