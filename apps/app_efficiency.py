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
import plotly.express as px
from flask_caching import Cache

import lib, header

from app_holder import app

# Import 
coef_fold = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/'
sys.path.append(coef_fold)
coef_df = pd.read_csv(coef_fold + f'assets/efficiency_coefficients.csv', index_col = 0)

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
        id = 'E_hidden_pred_signal_bystander',
        children = 'init'
      ),

      dcc.Location(
        id = 'E_url',
        refresh = False,
      ),

      html.Button(
        id = 'E_hidden_button_detect_pageload',
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
    header.get_navigation_header('efficiency'),

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
          html.Strong('Comparison of mean base editing efficiency by base editor')
          ],
          className = 'module_header_text'),
        ],
        className = 'module_header'
      ),

      # # Row
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
            transform = 'translate(35px, 115px)',
            fontSize = '14px',
          ),
        ),

        # Item. Plot
        dcc.Graph(
          id = 'E_efficiency_plot',
          config = dict(
            modeBarButtonsToRemove = modebarbuttons_2d,
            displaylogo = False,
            displayModeBar = False,
            staticPlot = True,
          ),
          style = dict(
            transform = 'translateX(20px)',
          ),
          className = 'ten columns',
        ),

        ],
        className = 'row',
        style = dict(
          marginBottom = '15px',
        ),
      ),

      # # Row: Slider
      html.Div([
        # Item. 
        html.Div(
          'Rescale editing efficiencies:',
          className = 'one column',
          style = dict(
            width = '32.4%',
            textAlign = 'right',
            transform = 'translate(-10px, 8px)',
          ),
        ),

        # Item
        html.Div(
          dcc.Slider(
            id = 'E_slider_efficiency_mean',
            min = 0.01,
            max = 0.99,
            step = 0.01,
            value = 0.50,
            updatemode = 'drag',
            # updatemode = 'mouseup',
          ),
          style = dict(
            float = 'left',
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

      # # Explanation text
      html.Div([
        # Item. 
        html.Div(
          [
            html.Div(
              'Mean base editing efficiencies were calculated in genome-integrated libraries of 12,000 sgRNA-target sites pairs designed with all 4-mers from protospacer positions 1-12. We performed logistic regression with a design matrix including covariates for 24 experimental batches, 3 cell-types, and each base editor. The fitted model has an R-squared value of 0.862. This page is provided to relate average base editing efficiencies across editors. Importantly, editing efficiency at any single target site cannot be converted across editors using this page alone - sequence preferences also differ across editors. To obtain a predicted efficiency at a single target site across editors, you should first convert mean editing efficiency between editors, then plug in the mean editing efficiency into the single mode efficiency module.'
            ),
          ],
          style = dict(
            transform = 'translate(120px, 0px)',
            marginBottom = '25px',
            width = '600px',
          ),
        ),

        ],
        className = 'row',
        style = dict(
        ),
      ),

    ], className = 'module_style',
    ),

    ],
    # body style
    # id = 'E_plots_body',
    style = dict(
      # display = 'none',
      transform = 'translateY(130px)',
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

###########################################
########     Module callbacks     #########
###########################################

##
# Efficiency
##
@app.callback(
  Output('E_efficiency_plot', 'figure'),
  [Input('E_slider_efficiency_mean', 'value'),
  ])
def efficiency_logit_plot(mean):
  from scipy.special import logit, expit
  logit_mean = logit(mean)
  std = 1
  # std = lib.efficiency_model_std

  # Form curve
  curve_xs = np.arange(-2.2, 2.2, 0.01)
  curve_ys = expit(std * curve_xs + logit_mean)

  all_data = [
    go.Scattergl(
      x = curve_xs,
      y = curve_ys,
      line = dict(
        color = lib.rgb['gray'],
      ),
      mode = 'lines',
      showlegend = False,
    ),
  ]

  special_pts_x, special_pts_y = [], []
  pt_colors = []
  pt_texts = []

  editor_shapes = []
  for idx, row in coef_df.iterrows():
    editor = row['Public base editor']
    x_val = float(row['Coefficient'])
    y_val = expit(x_val + logit_mean)
    color = lib.editor_cmap[editor]

    special_pts_x.append(x_val)
    special_pts_y.append(y_val)
    pt_colors.append(color)
    pt_texts.append(f'{editor}')
    # pt_texts.append(f'{editor}: {y_val:.1%}')

    editor_shapes.append(
      dict(
        type = 'line',
        x0 = x_val,
        y0 = 0,
        x1 = x_val,
        y1 = y_val,
        line = dict(
          color = color,
          width = 1.75,
          dash = 'dot',
        ),
      ),
    )

    all_data.append(
      go.Scatter(
        x = [x_val],
        y = [y_val],
        # mode = 'markers+text',
        mode = 'markers',
        marker = dict(
          color = color,
          size = 10,
        ),
        textposition = 'bottom center',
        hoverinfo = 'x+y',
        name = f'{editor}: {y_val:.0%}',
      ),
    )

  return dict(
    data = all_data, 
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
        title = 'Batch-adjusted logistic regression coefficient',
        tickvals = [-2, -1, 0, 1, 2],
        zeroline = False,
      ),
      shapes = editor_shapes, 
      height = 300,
      # width = 550,
      width = 700,
      margin = dict(
        l = 30,
        r = 30,
        t = 30,
        b = 30,
      ),
    ),
  )



