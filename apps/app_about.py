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

              This web app was developed by [Max W. Shen](https://www.maxwshen.com) using PyTorch, Dash Plotly, Heroku, and GitHub. Max W. Shen also developed the [inDelphi](https://www.crisprindelphi.design) online interactive web app which assists in designing CRISPR experiments. 

              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(
                  src = '/assets/graphical_abstract.png',
                  style = dict(
                    maxWidth = '70%',
                  ),
                ),
              ],
              style = dict(
                textAlign = 'center',
              )
            ),

            dcc.Markdown(dedent('''
              Copyright 2020. The Broad Institute, Massachusetts Institute of Technology, and President and Fellows of Harvard College. All rights reserved, except those expressly identified below.

              ***Limited License and Terms of Use for Research Use by Non-Profit and Government Institutions***
              
              As used in this limited license and terms of use, “Code” means this website, all software code running this website, all software code and models associated with this website, and all software code and models associated the following manuscript, no matter how the Code is accessed or obtained:

              __Mandana Arbab\*, Max W. Shen\*__, Beverly Mok, Christopher Wilson, Żaneta Matuszek, Christopher A. Cassa, and David R. Liu. "Determinants of Base Editing Outcomes from Target Library Analysis and Machine Learning." _Cell_, 2020, in press.

              BY DOWNLOADING OR USING THE CODE, YOU ARE CONSENTING TO BE AND AGREE TO BE BOUND BY ALL OF THE TERMS OF THIS LIMITED LICENSE, WHICH FOLLOW:

              • You agree you will not use the Code, or any information derived from the Code, for any commercial purposes whatsoever.

              • You agree you are an actively enrolled student, post-doctoral researcher, or faculty member at a degree-granting educational institution or US government research institution. 

              • You agree you will only use the Code for educational, instructional, and/or non-commercial research purposes. 

              • You understand that to obtain any right to use the Code for commercial purposes, including in the context of industrially sponsored research, you must enter into an appropriate separate and direct license agreement. For information on such licensing, please contact The Broad Institute using the following link: https://www.broadinstitute.org/contact 

              • You agree you will not redistribute unmodified versions of the Code.

              • You agree you will only modify the Code for educational, instructional, and/or non-commercial research purposes. 

              • You agree that if you redistribute any modifications of the Code, you will redistribute those modifications under the same terms as this license and only to non-profits and US government institutions.

              • For any use of the Code, including any modifications of the code, you agree you will credit the authors of the following manuscript: 

               __Mandana Arbab\*, Max W. Shen\*__, Beverly Mok, Christopher Wilson, Żaneta Matuszek, Christopher A. Cassa, and David R. Liu. "Determinants of Base Editing Outcomes from Target Library Analysis and Machine Learning." _Cell_, 2020, in press.

              • You agree to use neither the names of the owners (The Broad Institute, Massachusetts Institute of Technology (“MIT”), and President and Fellows of Harvard College (“Harvard”) or the names of the authors of the manuscript (above), to endorse or promote products or information derived from this software, without specific, prior written permission.

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


