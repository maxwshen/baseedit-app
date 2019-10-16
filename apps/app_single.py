import pickle, copy, os, datetime, subprocess, json, sys
from collections import defaultdict
import random
import numpy as np
import pandas as pd
import scipy
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
        id = 'S_hidden_pred_signal_bystander',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden_pred_signal_efficiency',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden_chosen_base_editor',
        children = 'BE4',
      ),
      html.Div(
        id = 'S_hidden_chosen_celltype',
        children = 'mES',
      ),
      html.Div(
        id = 'S_hidden_chosen_aa_frame',
        children = 'None',
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
  html.Div([
    ###################################################
    # Upper header
    ###################################################
    header.get_navigation_header('single'),

    ###################################################
    # Row: Sequence boxes
    ###################################################
    html.Div([
      # Item
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

    ###################################################
    # Row: Options
    ###################################################
    html.Div([
      # Item
      html.Div('',
        className = 'three columns',
      ),

      # Item
      dcc.Dropdown(
        id = 'S_editor_dropdown',
        options = lib.editor_celltype_dropdown_options,
        value = 'BE4, mES',
        searchable = True,
        clearable = False,
        style = dict(
          fontFamily = 'monospace',
          fontSize = 16,
        ),
        className = 'three columns'
      ),

      # Item
      dcc.Dropdown(
        id = 'S_aa_frame_dropdown',
        options = [
          {'label': 'None', 'value': 'None'},
          {'label': 'Frame 1, + strand', 'value': '1,+'},
          {'label': 'Frame 2, + strand', 'value': '2,+'},
          {'label': 'Frame 3, + strand', 'value': '3,+'},
          {'label': 'Frame 1, - strand', 'value': '1,-'},
          {'label': 'Frame 2, - strand', 'value': '2,-'},
          {'label': 'Frame 3, - strand', 'value': '3,-'},
        ],
        value = 'None',
        searchable = True,
        clearable = False,
        style = dict(
          fontFamily = 'monospace',
          fontSize = 16,
        ),
        className = 'three columns'
      ),

      # Item
      html.Div('',
        className = 'three columns',
      ),

      ], 
      style = dict(
        marginBottom = '5px',
        marginTop = '10px',
      ),
      className = 'row',
    ),

    # Row
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
  html.Div([
    ###################################################
    # Module: Efficiency
    ###################################################
    html.Div([
      # header
      html.Div([
        html.Div([
          html.Strong('Base editing efficiency')
          ],
          className = 'module_header_text'),
        ],
        className = 'module_header'
      ),

      # Row
      html.Div([
        # Item. Y-axis label
        html.Div(
          'Predicted frequency of sequenced reads with base editing activity at any substrate nucleotide',
          className = 'two columns',
          style = dict(
            textAlign = 'right',
            lineHeight = '1',
            # height = '200px',
            # translateY 100px sets top of text to middle, which makes multi-line text too low
            transform = 'translate(35px, 65px)',
            fontSize = '14px',
          ),
        ),

        # Item. Plot
        dcc.Graph(
          id = 'S_efficiency_plot',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
            staticPlot = True,
          ),
          style = dict(
            transform = 'translateX(20px)',
          ),
          className = 'four columns',
        ),

        # Item. Text descriptions
        html.Div(
          id = 'S_efficiency_longtext',
          style = dict(
            transform = 'translateY(15px)',
          ),
          className = 'six columns',
        ),

        ],
        className = 'row',
      ),

      # Row
      html.Div([
        # Item. 
        html.Div(
          '',
          className = 'three columns',
        ),

        # Item
        html.Div(
          dcc.Slider(
            id = 'S_slider_efficiency_mean',
            min = 0.01,
            max = 0.99,
            step = 0.01,
            value = 0.30,
            updatemode = 'drag',
            # updatemode = 'mouseup',
            marks = {
              0.10: {'label': '10%', 
                'style': {'color': lib.rgb['gray']},
              },
              0.25: {'label': '25%',
                'style': {'color': lib.rgb['gray']},
              },
              0.50: {'label': '50%',
                'style': {'color': lib.rgb['gray']},
              },
              0.75: {'label': '75%',
                'style': {'color': lib.rgb['gray']},
              },
              0.90: {'label': '90%',
                'style': {'color': lib.rgb['gray']},
              },
            },
          ),
          style = dict(
            float = 'right',
            width = '300px',
            marginTop = '15px',
            marginRight = '100px',
            marginBottom = '30px',
          ),
        ),

        ],
        className = 'row',
        style = dict(
          marginBottom = '10px',
        ),
      ),

    ], className = 'module_style',
    ),

    ###################################################
    # Module: Bystander, DNA
    ###################################################
    html.Div([
      # header
      html.Div([
        html.Div([
          html.Strong('Base editing outcomes: DNA sequence')
          ],
          className = 'module_header_text'),
        ],
        className = 'module_header'
      ),

      # Row
      html.Div([
        # Item
        # Text table
        dcc.Graph(
          id = 'S_bystander_gt_table',
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

    ###################################################
    # Module: Bystander, AA
    ###################################################
    html.Div([
      html.Div([
        # header
        html.Div([
          html.Div([
            html.Strong('Base editing outcomes: Amino acid sequence')
            ],
            className = 'module_header_text'),
          ],
          className = 'module_header'
        ),

        # Row
        html.Div([
          # Item
          # Text table
          dcc.Graph(
            id = 'S_bystander_aa_table',
            config = dict(
              modeBarButtonsToRemove = modebarbuttons_2d,
              displaylogo = False,
              displayModeBar = False,
              staticPlot = True,
            ),
            className = 'twelve columns',
          ),
          # Bar plot

          ],
          className = 'row',
        ),

        ], 
        className = 'module_style', 
      ),
      ], 
      id = 'S_bystander_module_container',
      style = {'display': 'none'},
      className = 'animate-bottom',
    ),

    ###################################################
    # Module
    ###################################################
    html.Div([
      ],
    )

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
# Hidden data callbacks
## 
@app.callback(
  Output('S_hidden_chosen_base_editor', 'children'),
  [Input('S_editor_dropdown', 'value')])
def update_editor_choice(val):
  [editor, celltype] = [s.strip() for s in val.split(',')]
  return editor

@app.callback(
  Output('S_hidden_chosen_celltype', 'children'),
  [Input('S_editor_dropdown', 'value')])
def update_celltype_choice(val):
  [editor, celltype] = [s.strip() for s in val.split(',')]
  return celltype

@app.callback(
  Output('S_hidden_chosen_aa_frame', 'children'),
  [Input('S_aa_frame_dropdown', 'value')])
def update_aaframe_choice(val):
  '''
    None
    1,+
    2,+
    3,+
    1,-
    2,-
    3,-
  '''
  return val

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
  Output('S_hidden_pred_signal_bystander', 'children'),
  [Input('S_textbox', 'value'),
   Input('S_hidden_chosen_base_editor', 'children'),
   Input('S_hidden_chosen_celltype', 'children')])
def bystander_predict(seq, base_editor, celltype):
  seq = seq.upper()
  bystander_predict_cache(seq, base_editor, celltype)
  return '%s,%s,%s' % (seq, base_editor, celltype)

@app.callback(
  Output('S_hidden_pred_signal_efficiency', 'children'),
  [Input('S_textbox', 'value'),
   Input('S_hidden_chosen_base_editor', 'children'),
   Input('S_hidden_chosen_celltype', 'children')])
def efficiency_predict(seq, base_editor, celltype):
  seq = seq.upper()
  efficiency_predict_cache(seq, base_editor, celltype)
  return '%s,%s,%s' % (seq, base_editor, celltype)

###########################################
########     Module callbacks     #########
###########################################

##
# Efficiency
##
@app.callback(
  Output('S_efficiency_longtext', 'children'),
  [Input('S_slider_efficiency_mean', 'value'),
   Input('S_hidden_pred_signal_efficiency', 'children'),
  ])
def update_efficiency_mean_text(chosen_mean, signal):

  seq, base_editor, celltype = signal.split(',')
  pred_d = efficiency_predict_cache(seq, base_editor, celltype)
  logit_score = pred_d['Predicted logit score']
  percentile = scipy.stats.norm.cdf(logit_score) * 100

  from scipy.special import logit, expit
  logit_mean = logit(chosen_mean)
  std = lib.efficiency_model_std
  pred_real = expit(std * logit_score + logit_mean) * 100

  if 0 <= percentile <= 35:
    var_text = 'below average'
    var_color = lib.rgb['red']
  if 35 <= percentile <= 65:
    var_text = 'average'
    var_color = 'black'
  if 65 <= percentile <= 100:
    var_text = 'above average'
    var_color = lib.rgb['green']

  # tooltip_msg = 'Base editing efficiency varies by cell-type, delivery method, length of exposure, etc. Here, we predict efficiency from sequence context alone.'
  tooltip_msg = 'Averaged across all possible sequence contexts'

  return [
    # Section
    html.Span('This target has '),
    html.Span(f'{var_text} base editing efficiency. ',
      style = dict(color = var_color),
    ),
    html.Br(),
    html.Span(f'Predicted Z-score: {logit_score:.2f}',
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span(f'Percentile: {percentile:.1f}',
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Br(),

    # Section
    html.Span(f'If the average editing efficiency is {100*chosen_mean:.1f}%, ',
    ),
    html.Div([
      html.Img(src = '/assets/tooltip_logo.png', className = 'tooltiplogo'),
      html.Span(tooltip_msg, className = 'tooltiptext'),
      ], className = 'tooltip',
    ),
    html.Br(),
    html.Span(f'then this target\'s efficiency is {pred_real:.1f}%.',
    ),
    html.Br(),
    html.Span(f'Adjust to your observed average editing efficiency in your experimental system. ',
      className = 'generalstats_subtext_style'),
    html.Span(f'Read more on how to do this accurately.',
      style = dict(fontStyle = 'italic'),
      className = 'generalstats_subtext_style'),
  ]

@app.callback(
  Output('S_efficiency_plot', 'figure'),
  [Input('S_slider_efficiency_mean', 'value'),
   Input('S_hidden_pred_signal_efficiency', 'children'),
  ])
def efficiency_logit_plot(mean, signal):
  seq, base_editor, celltype = signal.split(',')
  pred_d = efficiency_predict_cache(seq, base_editor, celltype)
  logit_score = pred_d['Predicted logit score']

  from scipy.special import logit, expit
  logit_mean = logit(mean)
  std = lib.efficiency_model_std

  # Form curve
  curve_xs = np.arange(-3, 3, 0.01)
  curve_ys = expit(std * curve_xs + logit_mean)

  # Form special points
  # special_pts_x = [0]
  # special_pts_y = [mean]
  special_pts_x = []
  special_pts_y = []
  special_pts_x.append(logit_score)
  pred_real = expit(std * logit_score + logit_mean)
  special_pts_y.append(pred_real)

  return dict(
    data = [
      go.Scattergl(
        x = curve_xs,
        y = curve_ys,
        line = dict(
          color = lib.rgb['blue'],
        ),
        mode = 'lines',
        hoverinfo = 'text',
        text = curve_ys,
        showlegend = False,
      ),
      go.Scattergl(
        x = special_pts_x,
        y = special_pts_y,
        mode = 'markers',
        marker = dict(
          color = lib.rgb['red'],
          size = 10,
        ),
        hoverinfo = 'x+y',
        showlegend = False,
      ),
    ],
    layout = go.Layout(
      font = dict(
        family = 'Arial',
      ),
      yaxis = dict(
        range = [0, 1],
        # title = 'Predicted fraction of sequenced reads with base editing activity at any substrate nucleotide',
        tickvals = [0, 0.25, 0.5, 0.75, 1],
        ticktext = ['0%', '25%', '50%', '75%', '100%'],
      ),
      xaxis = dict(
        title = 'Predicted efficiency score',
        tickvals = [-2, -1, 0, 1, 2],
        zeroline = False,
      ),
      shapes = [
        dict(
          type = 'line',
          x0 = logit_score,
          y0 = 0,
          x1 = logit_score,
          y1 = pred_real,
          line = dict(
            color = lib.rgb['red'],
            width = 1.75,
            dash = 'dot',
          ),
        ),
      ],
      height = 200,
      width = 300,
      margin = dict(
        l = 30,
        r = 30,
        t = 30,
        b = 30,
      ),
    ),
  )



##
# Bystander, genotype table
##
@app.callback(
  Output('S_bystander_gt_table', 'figure'),
  [Input('S_hidden_pred_signal_bystander', 'children'),
  ])
def update_gt_table(signal):
  seq, base_editor, celltype = signal.split(',')
  pred_df, stats, nt_cols = bystander_predict_cache(seq, base_editor, celltype)

  p0idx = 19
  target_seq = stats['50-nt target sequence']

  ## Set up data
  top10 = pred_df.iloc[:10]
  fqs = top10['Predicted frequency']
  fq_strings = [''] + [f'{100*s:.1f}%' for s in fqs]

  poswise_total = {col: sum(pred_df.loc[pred_df[col] != col[0], 'Predicted frequency']) for col in nt_cols}

  ## Form table with colors
  fillcolors = []
  fontcolors = []
  poswise_cols = []
  for gt_idx, ref_nt in enumerate(target_seq):
    pos = gt_idx - p0idx
    cand_col = f'{ref_nt}{pos}'
    pos_col = []
    col_fill_colors = []
    col_font_colors = []

    # row for target_seq 
    # Text
    pos_col.append(ref_nt)
    # Color
    if cand_col not in nt_cols:
      col_fill_colors.append('white')
    else:
      tot_edit_frac = poswise_total[cand_col]
      color_scale = lib.dna_color_scales[ref_nt]
      col_fill_colors.append(lib.get_color(color_scale, tot_edit_frac, white_threshold = 0.002))
    col_font_colors.append('black')

    # rows for edited genotypes
    for jdx, row in top10.iterrows():
      pred_fq = row['Predicted frequency']
      gt_seq = row['Genotype']

      obs_nt = gt_seq[gt_idx]

      # Genotype row
      # Text
      # Color
      pos_col.append(obs_nt)
      if obs_nt == ref_nt:
        col_fill_colors.append('white')
        col_font_colors.append(lib.font_cmap['match'])
      else:
        color_scale = lib.dna_color_scales[obs_nt]
        col_fill_colors.append(lib.get_color(color_scale, pred_fq))
        col_font_colors.append(lib.font_cmap['edited'])

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
# Bystander, amino acid table
##
@app.callback(
  Output('S_bystander_aa_table', 'figure'),
  [Input('S_hidden_pred_signal_bystander', 'children'),
   Input('S_hidden_chosen_aa_frame', 'children'),
  ])
def update_aa_table(signal, aa_frame_txt):
  seq, base_editor, celltype = signal.split(',')
  pred_df, stats, nt_cols = bystander_predict_cache(seq, base_editor, celltype)
  # if aa_frame_txt == 'None': return []
  aa_frame = int(aa_frame_txt[0]) - 1
  aa_strand = aa_frame_txt[-1]

  p0idx = 19
  target_seq = stats['50-nt target sequence']
  target_aas = lib.dna_to_aa(target_seq, aa_frame, aa_strand)
  pred_df['Amino acid sequence'] = [lib.dna_to_aa(s, aa_frame, aa_strand) for s in pred_df['Genotype']]
  aa_start_gap, aa_end_gap = lib.get_aa_display_start_idx(target_seq, aa_frame, aa_strand)

  # Whitespace ref AA
  whitespaced_target_aas = [' '] * aa_start_gap
  fills_target_aas = ['white'] * aa_start_gap
  for aa in target_aas:
    whitespaced_target_aas.append(f' {aa} ')
    fills_target_aas += [lib.aa_cmap[aa]] * 3
  whitespaced_target_aas += [' '] * aa_end_gap
  fills_target_aas += ['white'] * aa_end_gap
  whitespaced_target_aas = ''.join(whitespaced_target_aas)

  ## Set up data
  aa_fq_df = pred_df[['Amino acid sequence', 'Predicted frequency']].groupby('Amino acid sequence').agg(sum).reset_index().sort_values(by = 'Predicted frequency', ascending = False)
  aa_fq_df = aa_fq_df.iloc[:10]
  aa_fq_df = aa_fq_df[aa_fq_df['Predicted frequency'] >= 0.01]

  # Form amino acid whitespaced seqs, fill colors, font colors
  aa_to_fq = {}
  aa_to_gts = {}
  aa_to_fillcolors = {}
  aa_to_fontcolors = {}
  num_rows = 0
  for idx, row in aa_fq_df.iterrows():
    aa_seq = row['Amino acid sequence']

    whitespaced_aa_seq = [' '] * aa_start_gap
    fills_aas = ['white'] * aa_start_gap
    fontcolors_aas = ['white'] * aa_start_gap
    for jdx, aa in enumerate(aa_seq):
      target_aa = target_aas[jdx]
      if aa == target_aa:
        font_color = lib.font_cmap['match']
        fill_color = 'white'
      else:
        font_color = lib.font_cmap['edited']
        fill_color = lib.aa_cmap[aa]
      whitespaced_aa_seq.append(f' {aa} ')
      fills_aas += [fill_color] * 3
      fontcolors_aas += [font_color] * 3
    whitespaced_aa_seq += [' '] * aa_end_gap
    fills_aas += ['white'] * aa_end_gap
    fontcolors_aas += ['white'] * aa_end_gap
    whitespaced_aa_seq = ''.join(whitespaced_aa_seq)

    aa_to_fq[whitespaced_aa_seq] = row['Predicted frequency']
    num_rows += 2

    aa_to_fillcolors[whitespaced_aa_seq] = fills_aas
    aa_to_fontcolors[whitespaced_aa_seq] = fontcolors_aas

    dfs = pred_df[pred_df['Amino acid sequence'] == aa_seq]
    dfs = dfs[dfs['Predicted frequency'] >= 0.005].sort_values(by = 'Predicted frequency', ascending = False)
    aa_to_gts[whitespaced_aa_seq] = list(zip(list(dfs['Genotype']), list(dfs['Predicted frequency'])))
    num_rows += len(dfs)

  poswise_total = {col: sum(pred_df.loc[pred_df[col] != col[0], 'Predicted frequency']) for col in nt_cols}

  ## Form table with colors
  fillcolors = []
  fontcolors = []
  poswise_cols = []
  for gt_idx, ref_nt in enumerate(target_seq):
    pos = gt_idx - p0idx
    cand_col = f'{ref_nt}{pos}'
    pos_col = []
    col_fill_colors = []
    col_font_colors = []

    # Ref. amino acid seq row
    # Text and text
    ref_aa = whitespaced_target_aas[gt_idx]
    pos_col.append(ref_aa)
    col_fill_colors.append(fills_target_aas[gt_idx])
    col_font_colors.append('black')

    # Ref. target_seq row
    # Text
    pos_col.append(ref_nt)
    # Color
    if cand_col not in nt_cols:
      col_fill_colors.append('white')
    else:
      tot_edit_frac = poswise_total[cand_col]
      color_scale = lib.dna_color_scales[ref_nt]
      col_fill_colors.append(lib.get_color(color_scale, tot_edit_frac, white_threshold = 0.002))
    col_font_colors.append('black')

    # Blank row
    pos_col.append('')
    for col in col_fill_colors, col_font_colors:
      col.append('white')

    # rows for edited genotypes + aas
    for aa_seq in aa_to_fq:
      aa_fq = aa_to_fq[aa_seq]

      # Amino acid row
      # Text and color
      pos_col.append(aa_seq[gt_idx])
      col_fill_colors.append(aa_to_fillcolors[aa_seq][gt_idx])
      col_font_colors.append(aa_to_fontcolors[aa_seq][gt_idx])

      for (gt_seq, gt_fq) in aa_to_gts[aa_seq]:
        obs_nt = gt_seq[gt_idx]

        # Genotype row
        # Text
        # Color
        pos_col.append(obs_nt)
        if obs_nt == ref_nt:
          col_fill_colors.append('white')
          col_font_colors.append(lib.font_cmap['match'])
        else:
          color_scale = lib.dna_color_scales[obs_nt]
          col_fill_colors.append(lib.get_color(color_scale, gt_fq))
          col_font_colors.append(lib.font_cmap['edited'])

      # Blank row
      pos_col.append('')
      for col in col_fill_colors, col_font_colors:
        col.append('white')

    # Finished iterating over one column
    poswise_cols.append(pos_col)
    fillcolors.append(col_fill_colors)
    fontcolors.append(col_font_colors)

  # Get frequency strings
  fq_strings = ['', '', '']
  fq_string_fontcolors = ['white', 'white', 'white']
  for aa_seq in aa_to_fq:
    aa_fq = aa_to_fq[aa_seq]

    fq_strings.append(f'{100*aa_fq:.1f}%')
    fq_string_fontcolors.append('black')

    for (gt_seq, gt_fq) in aa_to_gts[aa_seq]:
      fq_strings.append(f'{100*gt_fq:.1f}%')
      fq_string_fontcolors.append(lib.rgb['gray'])

    # Blank row
    fq_strings.append('')
    fq_string_fontcolors.append('white')

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
          color = fontcolors + [fq_string_fontcolors],
        ),
        height = 20,
      ),
    )],
    layout = go.Layout(
      font = dict(
        family = 'monospace',
      ),
      height = 60 + 20 * num_rows + 20,
      width = 629,
      margin = dict(
        l = 10,
        r = 0,
        t = 5,
        b = 5,
      ),
    ),
  )

@app.callback(
  Output('S_bystander_module_container', 'style'),
  [Input('S_hidden_chosen_aa_frame', 'children')],
  [State('S_bystander_module_container', 'style')])
def show_hide_aa_module(aa_frame_text, prev_style):
  if aa_frame_text == 'None':
    prev_style['display'] = 'none'
  else:
    if 'display' in prev_style:
      del prev_style['display']
  return prev_style


##
# Download callbacks
##

##
# Flask serving
##

##
# Page link callback
##
