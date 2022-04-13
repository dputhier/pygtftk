import dash
import dash_bootstrap_components as dbc
import dash_defer_js_import as dji
import numpy as np
import pandas
import pandas as pd
import plotly.express as px
import plotly.io as pio
from dash import Input, Output, html, State
from dash import dcc

app = dash.Dash(external_stylesheets=[dbc.themes.YETI])

app.config.suppress_callback_exceptions = True

####################################################################################################
# Diagram template
####################################################################################################
# pio.templates.default = "ggplot2"
pio.templates.default = "simple_white"


####################################################################################################
# Function to load and prepare ologram table (barplot)
####################################################################################################

def loading_and_preparing_ologram_table_barplot(table_path):
    ####################################################################################################
    # loading table
    #####################################################################################################

    d = pd.read_csv(table_path,
                    sep="\t", header=0)

    d["feature_type"] = [x.replace(":", "\n") for x in d["feature_type"]]

    ####################################################################################################
    # Melting table
    #####################################################################################################

    dm = d.copy()

    # Create table for summed_bp_overlaps statistics
    #################################################
    data_ni_s = dm[['feature_type', 'summed_bp_overlaps_expectation_shuffled', 'summed_bp_overlaps_true']]
    maximum_s = data_ni_s[['summed_bp_overlaps_expectation_shuffled', 'summed_bp_overlaps_true']].max(axis=1)
    data_ni_s.columns = ['Feature', 'Shuffled', 'True']

    fc_s = data_ni_s['True'] / (data_ni_s['Shuffled'] + 1)
    fc_s = fc_s.to_list()

    dmm_s = data_ni_s.melt(id_vars='Feature')
    dmm_s = dmm_s.assign(Statistic=['Total overlap length per region type'] * dmm_s.shape[0])
    dmm_s.columns = ['Feature', 'Type', 'Value', 'Statistic']
    dmm_s = dmm_s.assign(Variance=np.sqrt(dm['summed_bp_overlaps_variance_shuffled']))

    # P-value
    ##########

    text_s = dm.loc[data_ni_s.index, 'summed_bp_overlaps_pvalue']

    # Format the text
    #################

    def format_p_value(x):
        if x == 0.0:
            r = 'p<1e-320'  # If the p-value is ~0 (precision limit), say so
        elif x == -1:
            r = 'p=NA'  # If the p-value was -1, we write 'Not applicable'
        else:
            r = '' + 'p={0:.2g}'.format(x)  # Add 'p=' before and format the p value
        return r

    text_s = [format_p_value(p) for p in text_s]
    Pval_text_s = pandas.DataFrame(text_s)

    dmm_s = dmm_s.assign(Pval_1=dm['summed_bp_overlaps_pvalue'])
    dmm_s = dmm_s.assign(Pval_2=Pval_text_s)
    dmm_s = dmm_s.assign(Neg_binom=dm['summed_bp_overlaps_negbinom_fit_quality'])

    # Create table for nb_intersections statistics
    ###############################################

    data_ni_n = dm[['feature_type', 'nb_intersections_expectation_shuffled', 'nb_intersections_true']]
    maximum_n = data_ni_n[['nb_intersections_expectation_shuffled', 'nb_intersections_true']].max(axis=1)
    data_ni_n.columns = ['Feature', 'Shuffled', 'True']

    fc_n = data_ni_n['True'] / (data_ni_n['Shuffled'] + 1)
    fc_n = fc_n.to_list()

    dmm_n = data_ni_n.melt(id_vars='Feature')
    dmm_n = dmm_n.assign(Statistic=['Total nb. of intersections per region type'] * dmm_n.shape[0])
    dmm_n.columns = ['Feature', 'Type', 'Value', 'Statistic']
    dmm_n = dmm_n.assign(Variance=np.sqrt(dm['nb_intersections_variance_shuffled']))

    # P-value
    ########

    text_n = dm.loc[data_ni_n.index, 'nb_intersections_pvalue']
    text_n = [format_p_value(p) for p in text_n]
    Pval_text_n = pandas.DataFrame(text_n)

    dmm_n = dmm_n.assign(Pval_1=dm['nb_intersections_pvalue'])
    dmm_n = dmm_n.assign(Pval_2=Pval_text_n)
    dmm_n = dmm_n.assign(Neg_binom=dm['nb_intersections_negbinom_fit_quality'])

    # Merge s and n tables
    ######################

    dmm = dmm_n.append(dmm_s)

    return dmm


# user_table_path = "/Users/puthier/Documents/git/project_dev/pygtftk/pygtftk/data/hg38_chr1/H3K36me3_ologram_stats.tsv"
user_table_path = "~/Documents/projet_bioinfo/pygtftk/ologram_output/00_ologram_stats.tsv"
dmm = loading_and_preparing_ologram_table_barplot(user_table_path)


####################################################################################################
# Function to load and prepare ologram table (volcano plot)
####################################################################################################

def loading_and_preparing_ologram_table_volcano(table_path):
    ####################################################################################################
    # loading table
    #####################################################################################################

    d = pd.read_csv(table_path,
                    sep="\t", header=0)

    d["feature_type"] = [x.replace(":", "\n") for x in d["feature_type"]]

    # Preparing a dataframe containing N statistics
    #################################################
    mat_n = d[['feature_type',
               'nb_intersections_log2_fold_change',
               'nb_intersections_pvalue']]

    # Unavailable p-value are discarded
    ####################################

    mat_n = mat_n.drop(mat_n[mat_n.nb_intersections_pvalue == -1].index)

    # Pval set to 0 are changed to  1e-320
    mat_n.loc[mat_n['nb_intersections_pvalue'] == 0, 'nb_intersections_pvalue'] = 1e-320
    mat_n = mat_n.assign(minus_log10_pvalue=list(-np.log10(list(mat_n.nb_intersections_pvalue))))
    mat_n.columns = ['Feature', 'log2(FC)', 'p-value', '-log10(pvalue)']
    mat_n = mat_n.assign(Statistic=['Total nb. of intersections per region type'] * mat_n.shape[0])

    # Preparing a dataframe containing S statistics
    #################################################

    mat_s = d[['feature_type',
               'summed_bp_overlaps_log2_fold_change',
               'summed_bp_overlaps_pvalue']]
    # Unavailable p-value are discarded

    mat_s = mat_s.drop(mat_s[mat_s.summed_bp_overlaps_pvalue == -1].index)
    # Pval set to 0 are changed to  1e-320
    mat_s.loc[mat_s['summed_bp_overlaps_pvalue'] == 0, 'summed_bp_overlaps_pvalue'] = 1e-320
    mat_s = mat_s.assign(minus_log10_pvalue=list(-np.log10(list(mat_s.summed_bp_overlaps_pvalue))))
    mat_s.columns = ['Feature', 'log2(FC)', 'p-value', '-log10(pvalue)']
    mat_s = mat_s.assign(Statistic=['Total overlap length per region type'] * mat_s.shape[0])

    ####################################################################################################
    # Preparing table for volcano
    #####################################################################################################

    df_volc = mat_n.append(mat_s)
    df_volc = df_volc.assign(Size=['3'] * df_volc.shape[0])

    return df_volc


# user_table_path = "/Users/puthier/Documents/git/project_dev/pygtftk/pygtftk/data/hg38_chr1/H3K36me3_ologram_stats.tsv"
user_table_path = "~/Documents/projet_bioinfo/pygtftk/ologram_output/00_ologram_stats.tsv"
df_volc = loading_and_preparing_ologram_table_volcano(user_table_path)

####################################################################################################
# App Layout
####################################################################################################

app.layout = dbc.Container(
    html.Div([
        dbc.Tabs(
            [
                dbc.Tab(label="Barplot", tab_id="tab-1"),
                dbc.Tab(label="Volcano Plot", tab_id="tab-2"),
                dbc.Tab(label="Table", tab_id="tab-3"),
                dbc.Tab(label="Original Table", tab_id="tab-4"),
            ],
            id="tabs",
            active_tab="tab-1",
        ),
        html.Div(id="content"),
    ]
    )
)

####################################################################################################
# Barplot Layout
####################################################################################################

available_statistics = dmm.Statistic.unique()
feature_type = dmm.Feature.unique()
results = dmm.Statistic

navbar_barplot = dbc.NavbarSimple(
    children=[
        dbc.Button(
            "Statistics",
            id="barplot-button-statistics",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Features",
            id="barplot-button-feature",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Font size",
            id="barplot-button-fontsize",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Themes",
            id="barplot-button-theme",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Thick Angle",
            id="barplot-button-tickangle",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Colors",
            id="barplot-button-color",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Orientation",
            id="barplot-button-orientation",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Ordering",
            id="barplot-button-ordering",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
    ],
    brand="Menu",
    brand_href="#",
    color="primary",
    dark=True,
)

barplot_statistics_menu = dcc.Dropdown(
    id="barplot_statistics_menu",
    options=available_statistics,
    value=available_statistics[0],
    clearable=False
)

barplot_feature_menu = dbc.Checklist(
    id="barplot_feature_menu",
    options=[{"label": x, "value": x} for x in feature_type],
    value=feature_type
)

barplot_fontsize_menu = dcc.Input(
    id="barplot_fontsize_menu",
    type='number',
    value=9,
    debounce=True,
    min=1,
    step=1)

barplot_theme_menu = dcc.Dropdown(
    id="barplot_theme_menu",
    options=list(pio.templates),
    value="simple_white",
    clearable=False)

barplot_tickangle_menu = dcc.Input(
    id="barplot_tickangle_menu",
    type='number',
    value=-45,
    debounce=True,
    min=-360,
    max=360,
    step=1)

barplot_colorshuffle_menu = dbc.Input(
    type="color",
    id="barplot_colorshuffle_menu",
    value="#FBBD04",
    style={"width": 75, "height": 50},
)

barplot_colortrue_menu = dbc.Input(
    type="color",
    id="barplot_colortrue_menu",
    value="#018CBA",
    style={"width": 75, "height": 50},
)

barplot_orientation_menu = dbc.RadioItems(
    options=[
        {"label": "Vertical", "value": 'v'},
        {"label": "Horizontal", "value": 'h'},
    ],
    value='h',
    id="barplot_orientation_menu",
)

barplot_ordering_menu = dcc.Dropdown(
    id="barplot_ordering_menu",
    options=[{"label": "p-value", "value": "Pval_1"},
             {"label": "Feature", "value": "Feature"}],
    value="Feature",
    clearable=False)

barplot_sorting_menu = dbc.RadioItems(
    options=[
        {"label": "Ascending", "value": True},
        {"label": "Decreasing", "value": False},
    ],
    value=True,
    id="barplot_sorting_menu",
    inline=True,
)

barplot_collapsed_menu = html.Div(
    id="barplot_collapsed_menu",
    children=[
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[barplot_feature_menu]))),
            id="barplot_feature_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[barplot_statistics_menu]))),
            id="barplot_statistics_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[barplot_fontsize_menu]))),
            id="barplot_fontsize_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[barplot_theme_menu]))),
            id="barplot_theme_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[barplot_tickangle_menu]))),
            id="barplot_tickangle_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(
                html.Div(children=[barplot_colorshuffle_menu, barplot_colortrue_menu]))),
            id="barplot_color_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[barplot_orientation_menu]))),
            id="barplot_orientation_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(
                html.Div(children=[barplot_ordering_menu, barplot_sorting_menu]))),
            id="barplot_ordering_menu_collapse",
            is_open=False,
        ),
    ])

barplot_layout = html.Div(children=[
    navbar_barplot,
    barplot_collapsed_menu,
    dcc.Graph(id="barplot"),

], id="html_div_barplot")


####################################################################################################
# Update functions for barplot
####################################################################################################

@app.callback(
    Output('barplot', "figure"),
    [Input("barplot_statistics_menu", "value"),
     Input("barplot_feature_menu", "value"),
     Input("barplot_fontsize_menu", "value"),
     Input("barplot_theme_menu", "value"),
     Input("barplot_tickangle_menu", "value"),
     Input("barplot_colortrue_menu", "value"),
     Input("barplot_colorshuffle_menu", "value"),
     Input("barplot_orientation_menu", "value"),
     Input("barplot_ordering_menu", "value"),
     Input("barplot_sorting_menu", "value")])
def update_graph(barplot_statistics_menu,
                 barplot_feature_menu,
                 bar_text_font_size,
                 barplot_theme_menu,
                 barplot_tickangle_menu,
                 barplot_colortrue_menu,
                 barplot_colorshuffle_menu,
                 barplot_orientation_menu,
                 barplot_ordering_menu,
                 barplot_sorting_menu):
    # Apply user-defined theming
    ################################
    pio.templates.default = barplot_theme_menu

    # Subset the dataset based on user selection
    # (N or S statistics)
    ############################################

    dmm_displayed = dmm[dmm["Statistic"] == barplot_statistics_menu]

    # Subset the dataset based on user selection
    # (the feature to display)
    ###########################################
    mask = dmm_displayed.Feature.isin(barplot_feature_menu)

    # Ensure "Shuffled" appear first.
    # This is important, later, to display the p-values
    ###################################################
    dmm_displayed = dmm_displayed.sort_values(by=["Type"])

    # Order based on used request.
    ###################################################
    dmm_displayed = dmm_displayed.sort_values(by=[barplot_ordering_menu], ascending=barplot_sorting_menu)

    # Prepare a bar diagram
    #######################

    orientation = barplot_orientation_menu

    if orientation == 'h':
        height = 50 * len(dmm_displayed[mask].Feature.unique())
        x = "Value"
        y = 'Feature'
        width = height / 3 + 400

    else:
        width = 80 * len(dmm_displayed[mask].Feature.unique())
        height = min(width / 3 + 400, 800)
        y = "Value"
        x = 'Feature'

    # error_y='Variance',
    fig = px.bar(dmm_displayed[mask],
                 x=x,
                 y=y,
                 orientation=orientation,
                 color="Type",
                 barmode="group",
                 title=barplot_statistics_menu,
                 height=height,
                 width=width,
                 color_discrete_map={
                     'Shuffled': barplot_colorshuffle_menu,
                     'True': barplot_colortrue_menu
                 })

    # Add p_value
    # Not that simple to ensure proper ordering...
    ###############################################
    # p-value to display
    pval_displayed = dmm_displayed[mask][dmm_displayed["Type"] == "Shuffled"].Pval_2.tolist()
    # x coordinates
    x_coord = range(len(pval_displayed))
    # y coordinates
    feature_ordering = dmm_displayed[mask][dmm_displayed["Type"] == "Shuffled"].Feature.tolist()
    y_val_shuffled = dict(zip(feature_ordering,
                              dmm_displayed[mask][dmm_displayed["Type"] == "Shuffled"].Value.tolist()))
    y_val_true = dict(zip(dmm_displayed[mask][dmm_displayed["Type"] == "True"].Feature.tolist(),
                          dmm_displayed[mask][dmm_displayed["Type"] == "True"].Value.tolist()))

    y_coord = [max(y_val_shuffled[x], y_val_true[x]) for x in feature_ordering]

    if orientation == 'h':
        tmp = x_coord
        x_coord = reversed(y_coord)
        y_coord = reversed(tmp)
        pval_displayed = reversed(pval_displayed)
        yshift = 0
        xshift = 30
    else:
        yshift = 20
        xshift = 0

    for x_val, y_val, label in zip(x_coord, y_coord, pval_displayed):
        fig.add_annotation(x=x_val, y=y_val, text=label, showarrow=False,
                           font={'size': int(bar_text_font_size), 'color': 'black'},
                           yshift=yshift,
                           xshift=xshift)

    fig.update_xaxes(tickangle=barplot_tickangle_menu)
    fig.update_layout(barmode='group')
    return fig


####################################################################################################
# Update functions for Menu / settings (barplot)
####################################################################################################

@app.callback(
    Output("barplot_feature_menu_collapse", "is_open"),
    [Input("barplot-button-feature", "n_clicks")],
    [State("barplot_feature_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_statistics_menu_collapse", "is_open"),
    [Input("barplot-button-statistics", "n_clicks")],
    [State("barplot_statistics_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_fontsize_menu_collapse", "is_open"),
    [Input("barplot-button-fontsize", "n_clicks")],
    [State("barplot_fontsize_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_theme_menu_collapse", "is_open"),
    [Input("barplot-button-theme", "n_clicks")],
    [State("barplot_theme_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_tickangle_menu_collapse", "is_open"),
    [Input("barplot-button-tickangle", "n_clicks")],
    [State("barplot_tickangle_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_color_menu_collapse", "is_open"),
    [Input("barplot-button-color", "n_clicks")],
    [State("barplot_color_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_orientation_menu_collapse", "is_open"),
    [Input("barplot-button-orientation", "n_clicks")],
    [State("barplot_orientation_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("barplot_ordering_menu_collapse", "is_open"),
    [Input("barplot-button-ordering", "n_clicks")],
    [State("barplot_ordering_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open

    return is_open


"""
## In progress
@app.callback(
    [Output("barplot_feature_menu_collapse", "is_open"),
     Output("barplot_statistics_menu_collapse", "is_open"),
     Output("barplot_fontsize_menu_collapse", "is_open"),
     Output("barplot_theme_menu_collapse", "is_open"),
     Output("barplot_tickangle_menu_collapse", "is_open"),
     Output("barplot_orientation_menu_collapse", "is_open"),
     Output("barplot_ordering_menu_collapse", "is_open"),
     ],
    [Input("barplot-button-feature", "n_clicks"),
     Input("barplot-button-statistics", "n_clicks"),
     Input("barplot-button-fontsize", "n_clicks"),
     Input("barplot-button-theme", "n_clicks"),
     Input("barplot-button-tickangle", "n_clicks"),
     Input("barplot-button-orientation", "n_clicks"),
     Input("barplot-button-ordering", "n_clicks"),
     ],
    [State("barplot_feature_menu_collapse", "is_open"),
     State("barplot_statistics_menu_collapse", "is_open"),
     State("barplot_fontsize_menu_collapse", "is_open"),
     State("barplot_theme_menu_collapse", "is_open"),
     State("barplot_tickangle_menu_collapse", "is_open"),
     State("barplot_orientation_menu_collapse", "is_open"),
     State("barplot_ordering_menu_collapse", "is_open")
     ]
)
def toggle_collapse_all_but_one(button_feature_n,
                                button_statistics_n,
                                button_fontsize_n,
                                button_theme_n,
                                button_tickangle_n,
                                button_orientation_n,
                                button_ordering_n,
                                button_feature_is_open,
                                button_statistics_is_open,
                                button_fontsize_is_open,
                                button_theme_is_open,
                                button_tickangle_is_open,
                                button_orientation_is_open,
                                button_ordering_is_open
                                ):
    if button_feature_n:
        print("button_feature_n pressed; state {button_feature_n}".format(button_feature_n=str(button_feature_n)))
        return True, False, False, False, False, False, False

    if button_statistics_n:
        print("button_feature_n pressed; state {button_feature_n}".format(button_feature_n=str(button_feature_n)))
        return False, True, False, False, False, False, False

    if button_statistics_n:
        print("button_feature_n pressed; state {button_feature_n}".format(button_feature_n=str(button_feature_n)))
        return False, True, False, False, False, False, False

"""

####################################################################################################
# Table tab
####################################################################################################

# The table layout
###################
table_layout = html.Div(children=[
    dbc.Card(
        dbc.CardBody([
            dbc.Table.from_dataframe(dmm,
                                     striped=True,
                                     bordered=True,
                                     hover=True,
                                     style={'margin-right': 'auto', 'margin-left': 'auto'},
                                     className="sortable"),
            dji.Import(src="https://www.kryogenix.org/code/browser/sorttable/sorttable.js")])
    )]
)

# The table layout
###################
original_table_layout = html.Div(children=[
    dbc.Card(
        dbc.CardBody([
            dbc.Table.from_dataframe(pd.read_csv(user_table_path,
                                                 sep="\t", header=0),
                                     striped=True,
                                     bordered=True,
                                     hover=True,
                                     style={'margin-right': 'auto', 'margin-left': 'auto'},
                                     className="sortable"),
            dji.Import(src="https://www.kryogenix.org/code/browser/sorttable/sorttable.js")])
    )]
)

####################################################################################################
# Volcano Layout
####################################################################################################
print(df_volc)
available_statistics_volc = df_volc.Statistic.unique()
feature_type_volc = df_volc.Feature.unique()
results_volc = df_volc.Statistic

navbar_volcano = dbc.NavbarSimple(
    children=[
        dbc.Button(
            "Statistics",
            id="volcano-button-statistics",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Features",
            id="volcano-button-feature",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Font size",
            id="volcano-button-fontsize",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Themes",
            id="volcano-button-theme",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Thick Angle",
            id="volcano-button-tickangle",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Colors",
            id="volcano-button-color",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
        dbc.Button(
            "Orientation",
            id="volcano-button-orientation",
            className="mb-1",
            color="primary",
            n_clicks=0,
        ),
    ],
    brand="Menu",
    brand_href="#",
    color="primary",
    dark=True,
)

volcano_statistics_menu = dbc.Checklist(
    id="volcano_statistics_menu",
    options=[{"label": x, "value": x} for x in available_statistics_volc],
    value=available_statistics_volc
)

volcano_feature_menu = dbc.Checklist(
    id="volcano_feature_menu",
    options=[{"label": x, "value": x} for x in feature_type_volc],
    value=feature_type_volc
)

volcano_fontsize_menu = dcc.Input(
    id="volcano_fontsize_menu",
    type='number',
    value=9,
    debounce=True,
    min=1,
    step=1)

volcano_theme_menu = dcc.Dropdown(
    id="volcano_theme_menu",
    options=list(pio.templates),
    value="simple_white",
    clearable=False)

volcano_tickangle_menu = dcc.Input(
    id="volcano_tickangle_menu",
    type='number',
    value=-45,
    debounce=True,
    min=-360,
    max=360,
    step=1)

volcano_N_menu = dbc.Input(
    type="color",
    id="volcano_N_menu",
    value="#FBBD04",
    style={"width": 75, "height": 50},
)

volcano_S_menu = dbc.Input(
    type="color",
    id="volcano_S_menu",
    value="#018CBA",
    style={"width": 75, "height": 50},
)

volcano_orientation_menu = dbc.RadioItems(
    options=[
        {"label": "Vertical", "value": 'v'},
        {"label": "Horizontal", "value": 'h'},
    ],
    value='h',
    id="volcano_orientation_menu",
)

volcano_collapsed_menu = html.Div(
    id="volcano_collapsed_menu",
    children=[
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[volcano_feature_menu]))),
            id="volcano_feature_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[volcano_statistics_menu]))),
            id="volcano_statistics_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[volcano_fontsize_menu]))),
            id="volcano_fontsize_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[volcano_theme_menu]))),
            id="volcano_theme_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[volcano_tickangle_menu]))),
            id="volcano_tickangle_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(
                html.Div(children=[volcano_N_menu, volcano_S_menu]))),
            id="volcano_color_menu_collapse",
            is_open=False,
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(html.Div(children=[volcano_orientation_menu]))),
            id="volcano_orientation_menu_collapse",
            is_open=False,
        ),
    ])

volcano_layout = html.Div(children=[
    navbar_volcano,
    volcano_collapsed_menu,
    dcc.Graph(id="volcano"),

], id="html_div_volcano")


####################################################################################################
# Update functions for volcano
####################################################################################################

@app.callback(
    Output('volcano', "figure"),
    [Input("volcano_statistics_menu", "value"),
     Input("volcano_feature_menu", "value"),
     Input("volcano_fontsize_menu", "value"),
     Input("volcano_theme_menu", "value"),
     Input("volcano_tickangle_menu", "value"),
     Input("volcano_N_menu", "value"),
     Input("volcano_S_menu", "value"),
     Input("volcano_orientation_menu", "value"),
     ])
def update_graph_volcano(volcano_statistics_menu,
                         volcano_feature_menu,
                         bar_text_font_size,
                         volcano_theme_menu,
                         volcano_tickangle_menu,
                         volcano_N_menu,
                         volcano_S_menu,
                         volcano_orientation_menu):
    # Apply user-defined theming
    ################################
    pio.templates.default = volcano_theme_menu

    # Subset the dataset based on user selection
    # (N or S statistics)
    ############################################

    df_volc_displayed = df_volc[df_volc.Statistic.isin(volcano_statistics_menu)]

    # Subset the dataset based on user selection
    # (the feature to display)
    ###########################################
    print("BLBAB")
    print(df_volc_displayed.Feature)
    print(volcano_feature_menu)
    mask_1 = df_volc_displayed.Feature.isin(volcano_feature_menu)
    print(mask_1)
    # Prepare a scatter plot diagram
    ################################

    orientation = volcano_orientation_menu

    if orientation == 'h':
        height = 500
        x = "Value"
        y = 'Feature'
        width = 500

    else:
        width = 500
        height = 500
        y = "Value"
        x = 'Feature'

    fig = px.scatter(df_volc_displayed[mask_1],
                     x='log2(FC)',
                     y='-log10(pvalue)',
                     color='Statistic',
                     hover_name='Feature')

    fig.update_traces(marker_coloraxis='coloraxis',
                      selector=dict(type='scatter'))

    x_coord = df_volc_displayed[mask_1]["log2(FC)"].tolist()
    y_coord = df_volc_displayed[mask_1]["-log10(pvalue)"].tolist()
    xy_label = df_volc_displayed[mask_1].Feature.tolist()

    for x_val, y_val, label in zip(x_coord, y_coord, xy_label):
        fig.add_annotation(x=x_val, y=y_val, text=label, showarrow=True)

    fig.update_xaxes(zeroline=True,
                     zerolinewidth=2,
                     zerolinecolor='LightPink')

    fig.update_yaxes(zeroline=True,
                     zerolinewidth=2,
                     zerolinecolor='LightPink')

    fig.update_traces(textposition='top center')

    fig.update_xaxes(tickangle=volcano_tickangle_menu)

    return fig


####################################################################################################
# Update functions for Menu / settings (Volcano)
####################################################################################################

@app.callback(
    Output("volcano_feature_menu_collapse", "is_open"),
    [Input("volcano-button-feature", "n_clicks")],
    [State("volcano_feature_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("volcano_statistics_menu_collapse", "is_open"),
    [Input("volcano-button-statistics", "n_clicks")],
    [State("volcano_statistics_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("volcano_fontsize_menu_collapse", "is_open"),
    [Input("volcano-button-fontsize", "n_clicks")],
    [State("volcano_fontsize_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("volcano_theme_menu_collapse", "is_open"),
    [Input("volcano-button-theme", "n_clicks")],
    [State("volcano_theme_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("volcano_tickangle_menu_collapse", "is_open"),
    [Input("volcano-button-tickangle", "n_clicks")],
    [State("volcano_tickangle_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("volcano_color_menu_collapse", "is_open"),
    [Input("volcano-button-color", "n_clicks")],
    [State("volcano_color_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    Output("volcano_orientation_menu_collapse", "is_open"),
    [Input("volcano-button-orientation", "n_clicks")],
    [State("volcano_orientation_menu_collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


####################################################################################################
# Update functions for Tabs
####################################################################################################

@app.callback(Output("content", "children"), [Input("tabs", "active_tab")])
def switch_tab(active_tab):
    ctx = dash.callback_context

    if active_tab == "tab-1":
        return barplot_layout
    elif active_tab == "tab-2":
        return volcano_layout
    elif active_tab == "tab-3":
        return table_layout
    elif active_tab == "tab-4":
        return original_table_layout
    else:
        return html.P("This shouldn't ever be displayed...")


app.run_server(debug=True, port=8071)
