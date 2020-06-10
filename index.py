import os

import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import flask

from app_holder import app
from apps import app_single
from apps import app_batch
from apps import app_efficiency
from apps import app_about
from apps import app_guide

base_dir = os.getcwd()

###################################################################
###################################################################
# Layout
app.layout = html.Div([
    html.Div(id = 'master-page-content'),
    dcc.Location(id = 'master-url', refresh = False),
])

app.title = 'BE-Hive'


###################################################################
###################################################################
# Serve pages
@app.callback(
  Output('master-page-content', 'children'),
  [Input('master-url', 'pathname')]
)
def display_page(pathname):
  # return app_single.layout
  print(pathname)
  if pathname is None or pathname == '/':
    return app_single.layout
  elif pathname[:len('/single')] == '/single':
    return app_single.layout
  elif pathname[:len('/batch')] == '/batch':
    return app_batch.layout
  elif pathname[:len('/efficiency')] == '/efficiency':
    return app_efficiency.layout
  elif pathname[:len('/guide')] == '/guide':
    return app_guide.layout
  elif pathname[:len('/about')] == '/about':
    return app_about.layout
  elif pathname[:len('/termsofuse')] == '/termsofuse':
    return app_termsofuse.layout
  else:
    return app_single.layout
  #   # return '404'

###################################################################
###################################################################

@app.server.route('/loading_dna.gif')
def flask_loading_gif():
  return flask.send_from_directory(base_dir + '/assets/', 'loading_dna.gif')

@app.server.route('/tooltip_logo.gif')
def flask_tooltip_logo():
  return flask.send_from_directory(base_dir + '/assets/', 'tooltip_logo.gif')

# Google analytics tracker
# Deprecated, use app_holder
# @app.server.route('/static/gtag.js')
# def serve_gtag():
#   return flask.send_from_directory(base_dir, 'gtag.js')

# app.scripts.append_script({'external_url': '/static/gtag.js'})

###################################################################
if __name__ == '__main__':
  app.run_server(debug = True)
  # app.run_server()