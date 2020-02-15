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
from collections import defaultdict

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
default_text = ''.join([random.choice(list('ACGT')) for s in range(100)])


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
        id = 'B_hidden_pred_signal_bystander',
        children = 'init'
      ),
      html.Div(
        id = 'B_hidden_pred_signal_efficiency',
        children = 'init'
      ),
      html.Div(
        id = 'B_hidden_chosen_base_editor_type',
        children = 'BE4',
      ),
      html.Div(
        id = 'B_hidden_chosen_celltype',
        children = 'mES',
      ),
      html.Div(
        id = 'B_hidden_chosen_aa_frame',
        children = 'None',
      ),

      dcc.Location(
        id = 'B_url',
        refresh = False,
      ),

      html.Button(
        id = 'B_hidden_button_detect_pageload',
        n_clicks_timestamp = 0,
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
    header.get_navigation_header('batch'),

    ###################################################
    # Row: Sequence boxes
    ###################################################
    html.Div([

      # Top row
      html.Div(
        [

          # Left half
          html.Div(
            [

              html.Strong(
                'Target genomic DNA: ',
                style = dict(
                  transform = 'translate(25%, 25%)',
                ),
                className = 'five columns',
              ),

              # Input text box: Target sequence
              dcc.Input(
                id = 'B_textbox', 
                size = '28',
                value = default_text,
                type = 'text',
                debounce = True,
                maxLength = 500,
                minLength = 50,
                persistence = True,
                autoFocus = True,
                style = dict(
                  fontFamily = 'monospace',
                  fontSize = 16,
                  # float = 'left',
                  # transform = 'translateX(-40px)',
                ),
                className = 'seven columns',
              ),
            ],
            className = 'six columns',
          ),

          # Right half
          html.Div(
            [
              html.Span(
                'Base editor / cell type:',
                style = dict(
                  transform = 'translate(35%, 25%)',
                ),
                className = 'six columns',
              ),

              # Item
              dcc.Dropdown(
                id = 'B_editor_dropdown',
                options = lib.batch_dropdown_options,
                value = lib.batch_dropdown_options[0]['value'],
                searchable = True,
                clearable = False,
                persistence = True,
                style = dict(
                  fontFamily = 'monospace',
                  fontSize = 16,
                ),
                className = 'six columns'
              ),
            ],
            className = 'six columns'
          ),

        ],
        className = 'row',
      ),

      # Bottom row
      html.Div(
        [
          # Left half
          html.Div(
            [
              html.Strong(
                'CRISPR protospacer: ',
                style = dict(
                  transform = 'translate(25%, 25%)',
                ),
                className = 'five columns',
              ),

              # Input text box: Protospacer
              dcc.Dropdown(
                id = 'B_protospacer_dropdown', 
                searchable = True,
                clearable = False,
                persistence = False,
                style = dict(
                  fontFamily = 'monospace',
                  fontSize = 16,
                  transform = 'translateX(10px)',
                  # transform = 'translateX(20px)',
                ),
                className = 'six columns',
              ),

              # Tooltip
              html.Div(
                [
                  html.Div(
                    [
                      html.Img(
                        src = 'https://www.crisprbehive.design/tooltip_logo.png',
                        className = 'tooltiplogo'
                      ),
                      html.Span(
                        "Must have ≥20 nt on 5' side and ≥10 nt on 3' side in target genomic DNA to be included", 
                        className = 'tooltiptext'
                      ),
                    ], 
                    style = dict(
                      transform = 'translate(10px, 8px)',
                    ),
                    className = 'tooltip',
                  ),
                ], 
                className = 'one columns'
              ),

            ],
            className = 'six columns'
          ),

          # Right half
          html.Div(
            [
              html.Span(
                'Amino acid frame:',
                style = dict(
                  transform = 'translate(35%, 25%)',
                ),
                className = 'six columns',
              ),

              dcc.Dropdown(
                id = 'B_aa_frame_dropdown',
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
                persistence = True,
                clearable = False,
                style = dict(
                  fontFamily = 'monospace',
                  fontSize = 16,
                ),
                className = 'five columns'
              ),

              html.Div('', 
                className = 'one columns'
              ),

            ],
            className = 'six columns',
          ),

        ],
        className = 'row',
      ),

      ], 
      style = dict(
        marginBottom = '5px',
        marginTop = '10px',
      ),
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
    # Module: Bystander, DNA
    ###################################################
    html.Div([
      # header
      html.Div([
        html.Div([
          html.Strong('Base editing outcomes among edited reads: DNA sequence')
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
          id = 'B_bystander_gt_table',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
            staticPlot = True,
          ),
          style = dict(
            height = 340,
            width = 800,
          ),
          className = 'twelve columns',
        ),

        ],
        className = 'row',
      ),

      # Download link
      html.Div(
        [
          html.Div(
            html.A(
              '📜 Download table of predictions', 
              id = 'B_csv_download_link'
            ),
          ),
        ], 
        style = dict(
          height = '30px',
          textAlign = 'center',
        ),
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
            html.Strong('Base editing outcomes among edited reads: Amino acid sequence')
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
            id = 'B_bystander_aa_table',
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
      id = 'B_bystander_module_container',
      style = {'display': 'none'},
      className = 'animate-bottom',
    ),

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
          id = 'B_efficiency_plot',
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
          id = 'B_efficiency_longtext',
          style = dict(
            transform = 'translateY(15px)',
          ),
          className = 'six columns',
        ),

        ],
        className = 'row',
      ),

      # Row: Slider
      html.Div([
        # Item. 
        html.Div(
          '',
          className = 'three columns',
        ),

        # Item
        html.Div(
          dcc.Slider(
            id = 'B_slider_efficiency_mean',
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

      html.Div(
        [
          html.Div(
            html.A(
              '🔗 Shareable link to your results', 
              id = 'B_page_link'
            ),
          ),
        ], 
        style = dict(
          height = '30px',
          textAlign = 'center',
        ),
      ),

    ], className = 'module_style',
    ),


    ###################################################
    # Module
    ###################################################
    html.Div([
      ],
    )

    ],
    # body style
    # id = 'B_plots_body',
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
# Input callbacks -- update with URL
##
@app.callback(
  Output('B_editor_dropdown', 'value'),
  [Input('B_hidden_button_detect_pageload', 'n_clicks_timestamp')],
  [State('B_url', 'pathname'),
   State('B_editor_dropdown', 'value')
  ],
)
def update_editor(timestamp, url, state):
  valid_flag, results_d = lib.parse_valid_url_path_single(url)
  if not valid_flag: 
    return state
  elif valid_flag:
    base_editor = results_d['base_editor']
    celltype = results_d['celltype']
    return f'{base_editor}, {celltype}'

@app.callback(
  Output('B_textbox', 'value'),
  [Input('B_hidden_button_detect_pageload', 'n_clicks_timestamp')],
  [State('B_url', 'pathname'),
   State('B_textbox', 'value')
  ],
)
def update_genomic_DNA(timestamp, url, state):
  valid_flag, results_d = lib.parse_valid_url_path_single(url)
  if not valid_flag: 
    return state
  elif valid_flag:
    return results_d['seq']

@app.callback(
  Output('B_protospacer_dropdown', 'value'),
  [Input('B_hidden_button_detect_pageload', 'n_clicks_timestamp')],
  [State('B_url', 'pathname'),
   State('B_textbox', 'value')
  ],
)
def update_protospacer(timestamp, url, seq):
  valid_flag, results_d = lib.parse_valid_url_path_single(url)
  if not valid_flag: 
    return f'{seq[20:40]}, 20'
  elif valid_flag:
    return results_d['protospacer']
    return f'{base_editor}, {celltype}'

@app.callback(
  Output('B_aa_frame_dropdown', 'value'),
  [Input('B_hidden_button_detect_pageload', 'n_clicks_timestamp')],
  [State('B_url', 'pathname'),
   State('B_aa_frame_dropdown', 'value')
  ],
)
def update_aa_frame(timestamp, url, state):
  valid_flag, results_d = lib.parse_valid_url_path_single(url)
  if not valid_flag: 
    return state
  elif valid_flag:
    return results_d['aa_frame_txt_parsed']
    return f'{base_editor}, {celltype}'

##
# Hidden data callbacks
## 
@app.callback(
  Output('B_hidden_chosen_base_editor_type', 'children'),
  [Input('B_editor_dropdown', 'value'),
  ])
def update_editor_choice(val):
  [editor_type, celltype] = [s.strip() for s in val.split(',')]
  return editor_type

@app.callback(
  Output('B_hidden_chosen_celltype', 'children'),
  [Input('B_editor_dropdown', 'value')])
def update_celltype_choice(val):
  [editor_type, celltype] = [s.strip() for s in val.split(',')]
  return celltype

@app.callback(
  Output('B_hidden_chosen_aa_frame', 'children'),
  [Input('B_aa_frame_dropdown', 'value')])
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
def bystander_predict_cache(seq, base_editor_type, celltype):
  bes = lib.type_to_bes[f'{base_editor_type}, {celltype}']
  mdf = None
  combined_stats = dict()
  for base_editor in bes:
    pred_df, stats = bystander_model.predict_given(
      seq,
      base_editor = base_editor,
      celltype = celltype,
    )
    pred_df = bystander_model.add_genotype_column(pred_df, stats)
    filtered_cols = ['Predicted frequency', 'Genotype']
    nt_cols = [col for col in pred_df if col not in filtered_cols]
    id_cols = nt_cols + ['Genotype']

    pred_df = pred_df.rename(columns = {
      'Predicted frequency': f'Predicted frequency, {base_editor}'
    })

    if mdf is None:
      mdf = pred_df
    else:
      mdf = mdf.merge(pred_df, on = id_cols, how = 'outer')

    # Combine stats (dicts)
    id_stats = [
      'Celltype',
      'Base editor',
    ]
    shared_stats = [
      '50-nt target sequence',
      'Assumed protospacer sequence',
    ]
    unique_stats = [
      'Total predicted probability',
    ]

    for ss in shared_stats:
      combined_stats[ss] = stats[ss]
    for us in unique_stats:
      combined_stats[f'{us}, {base_editor}, {celltype}'] = stats[us]

  # return pred_df, stats, nt_cols
  return mdf, combined_stats, nt_cols

##
# Prediction callbacks
##
@app.callback(
  Output('B_hidden_pred_signal_bystander', 'children'),
  [Input('B_textbox', 'value'),
   Input('B_protospacer_dropdown', 'value'),
   Input('B_hidden_chosen_base_editor_type', 'children'),
   Input('B_hidden_chosen_celltype', 'children')])
def bystander_predict(seq, ps, base_editor_type, celltype):
  seq = seq.upper()
  [protospacer, index] = ps.split()
  index = int(index)
  seq = seq[index - 20 : index + 30]
  bystander_predict_cache(seq, base_editor_type, celltype)
  return '%s,%s,%s' % (seq, base_editor_type, celltype)



##
# Protospacers
##
@app.callback(
  Output('B_protospacer_dropdown', 'options'),
  [Input('B_textbox', 'value')]
  )
def update_protospacers(seq):
  protospacers, poss, unique_flags = lib.find_protospacers(seq)

  options = []
  for protospacer, pos, unique_flag in zip(protospacers, poss, unique_flags):
    if unique_flag:
      d = {
        'label': f'{protospacer}',
        'value': f'{protospacer}, {pos}',
      }
    else:
      d = {
        'label': f'{protospacer}, {pos}',
        'value': f'{protospacer}, {pos}',
      }
    options.append(d)
  return options


###########################################
########     Module callbacks     #########
###########################################

##
# Bystander, genotype table
##
@app.callback(
  Output('B_bystander_gt_table', 'figure'),
  [Input('B_hidden_pred_signal_bystander', 'children'),
  ])
def update_gt_table(signal):
  seq, base_editor, celltype = signal.split(',')
  pred_df, stats, nt_cols = bystander_predict_cache(seq, base_editor, celltype)

  p0idx = 19
  target_seq = stats['50-nt target sequence']

  ## Set up data
  df_fq_cols = [col for col in pred_df.columns if 'Predicted frequency' in col]
  pred_df['Mean predicted frequency'] = pred_df[df_fq_cols].apply(np.nanmean, axis = 'columns')
  pred_df = pred_df.sort_values(by = 'Mean predicted frequency', ascending = False).reset_index(drop = True)
  top10 = pred_df.iloc[:10]

  fq_cols = []
  fq_fillcolors = []
  fq_fontcolors = []

  for df_fq_col in df_fq_cols:
    fq_strings = []
    curr_fills, curr_fonts = [], []

    editor_nm = df_fq_col.split(', ')[-1].strip()
    header = lib.editor_to_header[editor_nm]

    num_header_cols = 4
    header_strings = [''] * (num_header_cols - len(header))
    header_strings += header
    for _ in range(num_header_cols):
      fq_strings.append(header_strings[_])
      curr_fills.append('white')
      curr_fonts.append(lib.font_cmap['edited'])

    fqs = top10[df_fq_col]
    for fq in fqs:
      fq_strings.append(f'{100*fq:.0f}')
      if fq < 0.02:
        curr_fills.append('white')
        curr_fonts.append(lib.font_cmap['match'])
      else:
        curr_fills.append(lib.get_color(lib.gray_color_scale, fq))
        curr_fonts.append(lib.font_cmap['edited'])

    # fq_fillcolors.append(['white'] * len(fq_strings))
    # fq_fontcolors.append(['black'] * len(fq_strings))
    fq_cols.append(fq_strings)
    fq_fillcolors.append(curr_fills)
    fq_fontcolors.append(curr_fonts)


  poswise_total = {col: sum(pred_df.loc[pred_df[col] != col[0], 'Mean predicted frequency']) for col in nt_cols}

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

    # row for protospacer position numbers
    idx_to_num = {
      20: 1,
      24: 5,
      29: 1,
      30: 0,
      34: 1,
      35: 5,
      39: 2,
      40: 0,
    }
    if gt_idx in idx_to_num:
      pos_col.append(idx_to_num[gt_idx])
    else:
      pos_col.append('')
    col_font_colors.append('black')
    col_fill_colors.append('white')

    # row for |
    vert_idxs = [20, 24, 29, 34, 39]
    if gt_idx in vert_idxs:
      pos_col.append('|')
    else:
      pos_col.append('')
    col_font_colors.append('rgb(208, 211, 214)')
    col_fill_colors.append('white')

    # row for protospacer
    if gt_idx >= 20 and gt_idx < 40:
      pos_col.append(ref_nt)
    else:
      pos_col.append('')
    col_font_colors.append('black')
    col_fill_colors.append('white')

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
      pred_fq = row['Mean predicted frequency']
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
  fq_col_width = 5

  return dict(
    data = [go.Table(
      columnwidth = [pos_col_width] * len(poswise_cols) + [fq_col_width] * len(fq_cols),
      header = dict(
        line = dict(width = 0),
        fill = dict(color = 'white'),
        height = 0,
      ),
      cells = dict(
        values = poswise_cols + fq_cols,
        align = ['center'] * len(poswise_cols) + ['right'] * len(fq_cols), 
        fill = dict(
          color = fillcolors + fq_fillcolors,
        ),
        line = dict(width = 0),
        font = dict(
          family = 'monospace',
          color = fontcolors + fq_fontcolors,
        ),
        height = 20,
      ),
    )],
    layout = go.Layout(
      font = dict(
        family = 'monospace',
      ),
      width = 500 + 40 * len(fq_cols),
      margin = dict(
        l = 10,
        r = 0,
        t = 5,
        b = 5,
      ),
    ),
  )



