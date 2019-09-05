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
sys.path.append(app_fold + 'be_predict_bystander/')
import predict as bystander_model
sys.path[-1] = app_fold + 'be_predict_efficiency/'
import predict as efficiency_model

try:
  os.mkdir('user-csvs/')
except FileExistsError:
  pass
else:
  subprocess.check_output('rm -rf user-csvs/*', shell = True)

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

# Random default, which is cached on filesystem
default_left_text = ''.join([random.choice(list('ACGT')) for s in range(60)])
default_right_text = ''.join([random.choice(list('ACGT')) for s in range(4)]) + 'GG' + ''.join([random.choice(list('ACGT')) for s in range(54)])
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
        id = 'S_hidden-pred-signal',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden-cache-dsb-right',
        children = '%s' % (time.time())
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
  ),
  ##

  ],  # body div
  style = dict(
    width = '970px',
    margin = '0 auto',
  )
)


