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
    header.get_navigation_header('guide'),

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


    # Side bar
  html.Div(
    [
      ## 
      html.Div([
        html.A('Overview', 
          href = 'guide#overview', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(marginBottom = '8px')),

      html.Div([
        html.A('Python implementation', 
          href = 'guide#python', 
          className = 'dynamicunderline',
          style = dict(color = 'gray', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Calibrating editing efficiency predictions?', 
          href = 'guide#calibration', 
          className = 'dynamicunderline',
          style = dict(color = 'gray', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('What is total predicted probability?', 
          href = 'guide#tpp', 
          className = 'dynamicunderline',
          style = dict(color = 'gray', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Using BE-Hive with other Cas variants', 
          href = 'guide#casvariant', 
          className = 'dynamicunderline',
          style = dict(color = 'gray', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Using BE-Hive with other cell-types', 
          href = 'guide#celltype', 
          className = 'dynamicunderline',
          style = dict(color = 'gray', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

    ],
    style = dict(
      backgroundColor = 'white',
      opacity = 0.9,
      zIndex = 20,
      top = '12%',
      left = '2%',
      overflowX = 'hidden',
      position = 'fixed',
      width = '15%',
      # height = 300,
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
            dcc.Markdown(dedent('''
              BE-Hive is a suite of machine learning algorithms for assisting scientists using base editing. This interactive online web app is a companion to our publication:

              __Mandana Arbab\*, Max W. Shen\*__, Beverly Mok, Christopher Wilson, Żaneta Matuszek, Christopher A. Cassa, and David R. Liu. "Determinants of Base Editing Outcomes from Target Library Analysis and Machine Learning." _Cell_, 2020, in press.

              Please cite our paper if this web app was helpful in your work.

              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Overview

                  BE-Hive predicts base editing efficiency and bystander editing patterns for various base editors using machine learning models trained on observed base editing outcomes from up to 10,638 sgRNA-target sequence pairs integrated into the genomes of mouse embryonic stem cells and human HEK293T cells using SpCas9 and Cas9-NG base editors. These sgRNA-target pairs were designed to be minimally biased and maximally cover possible sequence space. Models for different base editors and cell-types were trained separately.

                  The input to any BE-Hive model is a genomic target sequence and an sgRNA sequence. The user selects which base editor and cell-type, which selects which machine learning models to use.

                  The editing efficiency model predicts the Z-score relative to the "average" sgRNA-target pair (across our dataset of highly diverse sgRNA-target pairs that cover sequence space with minimal bias). These Z-scores can be converted to the fraction of sequenced reads that have any base editing activity at any nucleotide in the base editing window among all sequenced reads, including unedited wild-type sequenced reads. (See calibration section below).

                  The bystander editing model predicts the frequency of a specific combination of base editing outcomes across all nucleotides in the base editing window among all sequenced reads that have any base editing activity at any nucleotide in the base editing window. 

                  Our single mode outputs predictions using the above units.

                  Predictions from the two models can be combined by simple multiplication since the units in the the bystander editing model's denominator and the editing efficiency model's numerator are the same. The units of the combined prediction are the frequency of a specific combination of base editing outcomes across all nucleotides in the base editing window among all sequenced reads, including unedited wild-type sequenced reads. Our batch mode combines predictions in this manner when the toggle "Report frequencies among: sequenced reads by including efficiency" is on.


                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'overview',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Python implementation

                  Python "packages" for the models are available at [https://github.com/maxwshen/be_predict_bystander](https://github.com/maxwshen/be_predict_bystander) and [https://github.com/maxwshen/be_predict_efficiency](https://github.com/maxwshen/be_predict_efficiency). If the capacity of the web application is insufficient for your needs, we recommend using these repos to run the BE-Hive models locally or on a server.


                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'python',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Calibrating editing efficiency predictions

                  Base editing efficiency depends on cell-type, delivery strategy, and other conditions unique to each experiment. To account for these factors, our base editing efficiency model outputs Z-scores by default, and allows users to provide experiment-specific information to convert the Z-score predictions to the units of the fraction of sequenced reads that have any base editing activity at any nucleotide in the base editing window among all sequenced reads, including unedited wild-type sequenced reads. 

                  The simplest strategy is to provide the "average" editing efficiency observed in your experimental system, where the average is taken over the theoretical set of all sgRNA-target pairs with all possible sequence contexts. Since most base editing experiments avoid sequence contexts known to have poor efficiency (such as those without centrally located cytosines when using cytosine base editors), simply averaging your previous base editing data is likely to overestimate this quantity.

                  In [https://github.com/maxwshen/be_predict_efficiency](https://github.com/maxwshen/be_predict_efficiency), we provide some simple python functions that will estimate this "average" quantity from observed data. 

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'calibration',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### What is total predicted probability?

                  In tables of predictions provided by our bystander editing model, there is a column called total predicted probability. 

                  The set of all possible combinations of editing outcomes grows exponentially with the number of substrate nucleotides in a base editing window (denote this as N). At a first glance, this number may appear to be 2^N when considering only two possibilities: that each single nucleotide is either edited or not (C or T, in the case of cytosine base editors). However, our work identifies uncommon and rare base editing outcomes including C->G, C->A, G->A conversions by CBEs. Thus when considering cytosine base editing, the possiblity space scales as 4^N. When N is large, it can take a prohibitive number of forward model evaluations to predict the probability of all exponentially many editing combinations, which would sum to 1.

                  However, we know that some editing combinations are more likely than others. To provide predictions in an efficient and expedient manner, we use greedy heuristics to minimize the number of forward model evaluations while maximizing the total probability accounted for. Since we only query the model on a subset of all possible sequences, the total probability observed must be less than 1. 

                  In typical cases, the total predicted probability is 0.95 or greater. For downstream applications, we recommend the conservative assumption that the remaining probability are allocated to the least desirable editing outcome possible.

                  We expect that in unusual cases, the total predicted probability provided using the default web app settings are insufficient for your needs. For example, a sequence context edited by cytosine base editors with many cytosines could have insufficiently high total predicted probability for your experimental design needs. Our python implementations at [https://github.com/maxwshen/be_predict_efficiency](https://github.com/maxwshen/be_predict_efficiency) allow finer-grained control over the total predicted probability, such as increasing runtime to increase total predicted probability.

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'tpp',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### How do I use BE-Hive with other cell-types?

                  We anticipate that base editing activity is generally similar across mammalian cell-types that share similar DNA repair systems. We recommend selecting between mES and HEK293T models by similarity of DNA repair systems. Another less important criteria that could be useful is that our mES data was generally cleaner and higher quality than our HEK293T data.

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'celltype',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### How do I use BE-Hive with other Cas variants?

                  Our web app does not explicitly filter protospacers by PAM. If your Cas variant has similar base editing activity as SpCas9 or Cas9-NG base editors, but has a different PAM, you can simply select the appropriate protospacers from the drop-down menus in our web app.

                  If your Cas variant base editor has different activity than SpCas9 or Cas9-NG base editors, including SaCas9 and Cas12a (Cpf1), please refer to our manuscript and supplementary information which discuss using BE-Hive trained on SpCas9/Cas9-NG base editing data on these Cas variants. Long story short, the base editing window tends to shift and sometimes widen or narrow when modifying the Cas variant, but deaminase-specific sequence preferences do not change substantially (as one would expect).

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'casvariant',
              className = 'hashbuffer',
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


