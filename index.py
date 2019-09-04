import os

import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import flask

from app import app
from apps import app_single, app_batch, app_gene, app_guide, app_about, app_termsofuse

###################################################################
###################################################################
# Layout
app.layout = html.Div([
    html.Div(id = 'master-page-content'),
    dcc.Location(id = 'master-url', refresh = False),

    # Hidden datatable for proper rendering
    # https://community.plot.ly/t/display-tables-in-dash/4707/40?u=chriddyp
    html.Div(dt.DataTable(rows=[{}]), style={'display': 'none'})
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
  elif pathname[:len('/gene')] == '/gene':
    return app_gene.layout
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
# CSS
css_directory = os.getcwd()
@app.server.route('/static/<stylesheet>')
def serve_stylesheet(stylesheet):
  if stylesheet not in stylesheets:
    raise Exception(
      '"{}" is excluded from the allowed static files'.format(
        stylesheet
      )
    )
  return flask.send_from_directory(css_directory, stylesheet)

stylesheets = ['stylesheet.css']
for stylesheet in stylesheets:
  app.css.append_css({'external_url': '/static/{}'.format(stylesheet)})

# As of 0.22.0, served automatically from /assets/

# Favicon
@app.server.route('/favicon.ico')
def favicon():
  return flask.send_from_directory(os.getcwd() + '/assets/', 'favicon.ico', mimetype = 'image/vnd.microsoft.icon')

# As of 0.22.0, served automatically from /assets/

# Google analytics tracker
@app.server.route('/static/gtag.js')
def serve_gtag():
  return flask.send_from_directory(css_directory, 'gtag.js')

app.scripts.append_script({'external_url': '/static/gtag.js'})

###################################################################
if __name__ == '__main__':
  app.run_server(debug = True)
  # app.run_server()