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
    ],
  ),

  ##
  # Header
  ##
  html.Div([
    ###################################################
    # Upper header
    ###################################################
    header.get_navigation_header('about'),

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
    ###################################################
    # Section: Overview
    ###################################################
    html.Div(
      [
        html.Div(
          [
            # html.Div(
            #   [
            #     html.Img(src = '/assets/fig-indel-len.gif'),
            #     html.Div(
            #       [
            #         ''
            #       ],
            #       style = dict(
            #         fontStyle = 'italic',
            #         width = '450px',
            #         margin = '0 auto',
            #         marginBottom = '20px',
            #       )
            #     ),
            #   ],
            #   style = dict(
            #     textAlign = 'center',
            #   ),
            # ),

            dcc.Markdown(dedent('''
              BE-Hive is a suite of machine learning algorithms for assisting scientists using base editing. 

              This interactive online web app is a companion to our publication:

              __Mandana Arbab\*, Max W. Shen\*__, Beverly Mok, Christopher Wilson, Żaneta Matuszek, Christopher A. Cassa, and David R. Liu. "Determinants of Base Editing Outcomes from Target Library Analysis and Machine Learning." _Cell_, 2020, in press.

              Please cite our paper if this web app was helpful in your work.

              This web app was developed by Max W. Shen using PyTorch, Dash Plotly, Heroku, and GitHub. Max W. Shen also developed the [inDelphi](www.crisprindelphi.design) online interactive web app which assists in designing CRISPR experiments. 

              '''),
              className = 'markdown_style',
            ),

            # html.Div(
            #   [
            #     html.Img(src = '/assets/fig-coverplus.PNG'),
            #     html.Div(
            #       [
            #         ''
            #       ],
            #       style = dict(
            #         fontStyle = 'italic',
            #         width = '450px',
            #         margin = '0 auto',
            #         marginBottom = '20px',
            #       )
            #     ),
            #   ],
            #   style = dict(
            #     textAlign = 'center',
            #   ),
            # ),

          ],
          style = dict(
            width = '800px',
            margin = '0 auto',
          ),
        ),
      ],
      id = 'overview',
      className = 'hashbuffer',
    ),
    ],
    style = dict(
      transform = 'translateY(120px)',
      marginBottom = '150px',
    )
  )
  ],  # body div
  style = dict(
    width = '1150px',
    margin = '0 auto',
  )
)


