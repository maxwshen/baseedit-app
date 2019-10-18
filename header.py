import dash
import dash_core_components as dcc
import dash_html_components as html
from app_holder import app

divider_text = '   •   '

def get_navigation_header(page_nm):
  font_size_param = 16

  dot_style = dict(
    color = 'gray',
    fontSize = '%spx' % (font_size_param),
  )

  default_style = dict(
    position = 'relative',
    textDecoration = 'none',
    textTransform = 'uppercase',
    fontFamily = 'sans-serif',
    color = 'black',
    marginBottom = '2px',
    fontSize = '%spx' % (font_size_param),
  )

  import copy
  selected_style = copy.copy(default_style)
  selected_style['borderBottom'] = '1px solid'

  styles = dict()
  for nm in ['single', 'batch', 'gene', 'guide', 'about']:
    if page_nm == nm:
      styles[nm] = selected_style
    else:
      styles[nm] = default_style

  return html.Div(
    [
      # html.H4(
      #   'BE-Hive',
      #   style = dict(
      #     textAlign = 'center',
      #   ),
      # ),
      html.Div(
        [
          html.Img(
            src = app.get_asset_url('bee.png'), 
            style = dict(
              height = '32px',
              verticalAlign = 'top',
              marginRight = 5,
            ),
          ),
          html.Span(
            'BE-Hive',
            style = dict(
              textAlign = 'center',
              fontSize = 26,
            ),
          ),
        ],
        style = dict(
          marginTop = 8,
          marginBottom = 8,
          textAlign = 'center',
        ),
      ),

      html.Div(
        [
          html.A(
            'Single mode',
            href = 'single',
            style = styles['single'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'Batch mode',
            href = 'batch',
            style = styles['batch'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'Gene mode',
            href = 'gene',
            style = styles['gene'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'User guide',
            href = 'guide',
            style = styles['guide'],
            className = 'dynamicunderline',
          ),
          html.Span(
            divider_text,
          ),
          html.A(
            'About',
            href = 'about',
            style = styles['about'],
            className = 'dynamicunderline',
          ),
        ],
        style = dict(
          marginBottom = 20,
          textAlign = 'center',
        ),
        className = 'row',
      ),

    ],
  )