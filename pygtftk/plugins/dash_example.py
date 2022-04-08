import dash
import numpy as np
import pandas
import pandas as pd
import plotly.express as px
import plotly.io as pio
from dash import dcc
from dash import html
from dash.dependencies import Input, Output

from pygtftk.utils import sort_2_lists

####################################################################################################
# Setting variables
#####################################################################################################

sort_features = "summed_bp_overlaps_true"
display_fit_quality = False
coord_flip = False
feature_order = None
only_those_combis = None

# Diagram template
###################
pio.templates.default = "ggplot2"

####################################################################################################
# loading table
#####################################################################################################

d = pd.read_csv("/Users/puthier/Documents/git/project_dev/pygtftk/pygtftk/data/hg38_chr1/H3K36me3_ologram_stats.tsv",
                sep="\t", header=0)

d["feature_type"] = [x.replace(":", "\n") for x in d["feature_type"]]

if sort_features is not None:
    sorted_feat = sort_2_lists(d[sort_features].tolist(),
                               d.feature_type.tolist())[1]
    feature_order = []
    for x in sorted_feat:
        if x not in feature_order:
            feature_order += [x]
else:
    feature_order = None

####################################################################################################
# Melting table
#####################################################################################################

dm = d.copy()

if only_those_combis is not None:
    dms = dm.loc[dm['feature_type'].isin(only_those_combis)]
else:
    dms = dm

# Create table for summed_bp_overlaps statistics
#################################################
data_ni_s = dms[['feature_type', 'summed_bp_overlaps_expectation_shuffled', 'summed_bp_overlaps_true']]
maximum_s = data_ni_s[['summed_bp_overlaps_expectation_shuffled', 'summed_bp_overlaps_true']].max(axis=1)
data_ni_s.columns = ['Feature', 'Shuffled', 'True']

fc_s = data_ni_s['True'] / (data_ni_s['Shuffled'] + 1)
fc_s = fc_s.to_list()

dmm_s = data_ni_s.melt(id_vars='Feature')
dmm_s = dmm_s.assign(Statistic=['Total overlap length per region type'] * dmm_s.shape[0])
dmm_s.columns = ['Feature', 'Type', 'Value', 'Statistic']
dmm_s = dmm_s.assign(Variance=np.sqrt(dms['summed_bp_overlaps_variance_shuffled']))

# P-value
##########

text_s = dms.loc[data_ni_s.index, 'summed_bp_overlaps_pvalue']


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

dmm_s = dmm_s.assign(Pval_1=dms['summed_bp_overlaps_pvalue'])
dmm_s = dmm_s.assign(Pval_2=Pval_text_s)
dmm_s = dmm_s.assign(Neg_binom=dms['summed_bp_overlaps_negbinom_fit_quality'])

# print(dmm_s)

# Create table for nb_intersections statistics
###############################################

data_ni_n = dms[['feature_type', 'nb_intersections_expectation_shuffled', 'nb_intersections_true']]
maximum_n = data_ni_n[['nb_intersections_expectation_shuffled', 'nb_intersections_true']].max(axis=1)
data_ni_n.columns = ['Feature', 'Shuffled', 'True']

fc_n = data_ni_n['True'] / (data_ni_n['Shuffled'] + 1)
fc_n = fc_n.to_list()

dmm_n = data_ni_n.melt(id_vars='Feature')
dmm_n = dmm_n.assign(Statistic=['Total nb. of intersections per region type'] * dmm_n.shape[0])
dmm_n.columns = ['Feature', 'Type', 'Value', 'Statistic']
dmm_n = dmm_n.assign(Variance=np.sqrt(dms['nb_intersections_variance_shuffled']))

# P-value
########

text_n = dms.loc[data_ni_n.index, 'nb_intersections_pvalue']
text_n = [format_p_value(p) for p in text_n]
Pval_text_n = pandas.DataFrame(text_n)

dmm_n = dmm_n.assign(Pval_1=dms['nb_intersections_pvalue'])
dmm_n = dmm_n.assign(Pval_2=Pval_text_n)
dmm_n = dmm_n.assign(Neg_binom=dms['nb_intersections_negbinom_fit_quality'])

# print(dmm_n)


# Merge s and n tables
######################

dmm = dmm_n.append(dmm_s)
dmm['Pval_2'] = dmm['Pval_2'].astype('string')
# print(dmm.dtypes)
# Order the categories in order to the p-value
###############################################

dmm = dmm.sort_values(by='Pval_1')
# dmm_2 = dmm
# print(dmm)
# print(dmm.dtypes)

# print(dmm.columns.tolist())

####################################################################################################
# Preparing table for barplot
#####################################################################################################

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

# print(df_volc)

# name = df_volc['Feature']
# name = str(name)
# print(name)

# x_scatter = df_volc['log2(FC)']
# x_scatter = str(x_scatter)
# print(x_scatter)

# y_scatter = df_volc['-log10(pvalue)']
# y_scatter = str(y_scatter)
# print(y_scatter)

####################################################################################################
# Plotting
#####################################################################################################

# Graph creation
#################

results_stat = dmm.Statistic.unique()
feature_type = dmm.Feature.unique()
results = dmm.Statistic
volcano = df_volc

app = dash.Dash(__name__)

volcano_plot = px.scatter(volcano, x='log2(FC)',
                          y='-log10(pvalue)',
                          color='Statistic',
                          hover_name='Feature',
                          text='Feature')

volcano_plot.update_traces(marker_coloraxis='coloraxis',
                           selector=dict(type='scatter'))

volcano_plot.update_xaxes(zeroline=True,
                          zerolinewidth=2,
                          zerolinecolor='LightPink')

volcano_plot.update_yaxes(zeroline=True,
                          zerolinewidth=2,
                          zerolinecolor='LightPink')

volcano_plot.update_traces(textposition='top center')

# volcano_plot.add_annotation(x='log2(FC)', y='-log10(pvalue)', text=name, showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2)


app.layout = html.Div([
    dcc.Dropdown(
        id="statistics_dropdown",
        options=results_stat,
        value=results_stat[0],
        clearable=False, ),
    dcc.Checklist(
        id="Feature_checklist",
        options=[{"label": x, "value": x} for x in feature_type],
        value=feature_type[0:], ),
    dcc.Input(
        id="bar_text_font_size",
        placeholder='P-value font size',
        type='text',
        value='14', ),
    #    dcc.Graph(id="chek"),
    dcc.Graph(id="barplot"),
    html.Hr(),
    dcc.Graph(id="scatter-plot", figure=volcano_plot),
    html.P("Legend position"),
    dcc.RadioItems(
        id='xanchor',
        options=[{'label': 'left', 'value': 0},
                 {'label': 'right', 'value': 1}],
        value=0,
        labelStyle={'display': 'inline-block'}
    ),
    dcc.RadioItems(
        id='yanchor',
        options=[{'label': 'top', 'value': 1},
                 {'label': 'bottom', 'value': 0}],
        value=1,
        labelStyle={'display': 'inline-block'}
    ),

])


# @app.callback(
#    Output("bar-chart", "figure"),
#    Input("dropdown", "value"),
#    Input("checklist", "value"))
# def update_bar_chart(Statistic, feature_list):
#    mask = dmm["Statistic"] == Statistic
#    fig = px.bar(dmm[mask], x="Feature", y="Value",
#                 color="Type", error_y='Variance', barmode="group", text='Pval_1')
#    #fig.add_trace(go.Bar(name='Shuffled', x='Feature', y='Value', error_y=dict(type='data', value='Variance)')))
#    fig.update_traces(texttemplate='%{text:.2s}', textposition='outside')
#    return fig


##################################
# if len(checklist) == 0:
#    chosen_feature = dmm_2[dmm_2['Feature'].isin(
#        ['Feature', 'Promoters', 'start_codon', 'five_prime_utr', 'exon', 'CDS', 'Introns', 'transcript', 'gene',
#         'Intergenic', 'Terminator', 'stop_codon'])]
# else:
#    chosen_feature = dmm_2[dmm_2.index.isin([checklist])]
# list_feature = list(chosen_feature['Feature'])

# print(list_feature)
# feat = dmm.Feature.isin(list_feature)
# mask = (feat["Statistic"] == Statistic)
# fig = px.bar(feat[mask], x="Feature", y="Value",
#             color="Type", error_y='Variance', barmode="group", text='Pval_1')
# fig.add_trace(go.Bar(name='Shuffled', x='Feature', y='Value', error_y=dict(type='data', value='Variance)')))
# fig.update_traces(texttemplate='%{text:.2s}', textposition='outside')
# return fig
###################################


# @app.callback(
#    Output("chek", "figure"),
#    Input("checklist", "value"))
# def update_chek(checklist):
#    mask = dmm.Feature.isin(checklist)
#    plot = px.bar(dmm[mask], x='Feature', y="Value",
#                  color="Type", error_y='Variance', barmode="group", text='Pval_1')
#    return plot

####################################################################################################
# Dynamic plots
#####################################################################################################

@app.callback(
    Output('barplot', "figure"),
    [Input("statistics_dropdown", "value"),
     Input("Feature_checklist", "value"),
     Input("bar_text_font_size", "value")])
def update_graph(statistics_dropdown, checklist, bar_text_font_size):
    # font size should be an int
    #############################
    if bar_text_font_size == '':
        bar_text_font_size = 14
    bar_text_font_size = int(bar_text_font_size)

    # Subset the dataset based on user selection
    # (N or S statistics)
    ############################################
    print("Selected Stat:" + statistics_dropdown)

    dmm_displayed = dmm[dmm["Statistic"] == statistics_dropdown]

    # Subset the dataset based on user selection
    # (the feature to display)
    ###########################################
    mask = dmm_displayed.Feature.isin(checklist)

    # Ensure "Shuffled" appear first.
    # This is important, later, to display the p-values
    ###################################################
    dmm_displayed = dmm_displayed.sort_values(by=["Type"])

    # Prepare a bar diagram
    #######################
    fig = px.bar(dmm_displayed[mask],
                 x='Feature',
                 y="Value",
                 color="Type",
                 error_y='Variance',
                 barmode="group",
                 title=statistics_dropdown)

    # Add p_value
    # Not that simple to ensure proper ordering...
    ###############################################
    # p-value to display
    pval_displayed = dmm_displayed[mask][dmm_displayed["Type"] == "Shuffled"]["Pval_2"].tolist()
    # x coordinates
    x_coord = range(len(pval_displayed))
    # y coordinates
    feature_ordering = dmm_displayed[mask][dmm_displayed["Type"] == "Shuffled"]["Feature"].tolist()
    y_val_shuffled = dict(zip(feature_ordering,
                              dmm_displayed[mask][dmm_displayed["Type"] == "Shuffled"]["Value"].tolist()))
    y_val_true = dict(zip(dmm_displayed[mask][dmm_displayed["Type"] == "True"]["Feature"].tolist(),
                          dmm_displayed[mask][dmm_displayed["Type"] == "True"]["Value"].tolist()))

    y_coord = [max(y_val_shuffled[x], y_val_true[x]) for x in feature_ordering]

    for x_val, y_val, label in zip(x_coord, y_coord, pval_displayed):
        fig.add_annotation(x=x_val, y=y_val, text=label, showarrow=False,
                           font={'size': int(bar_text_font_size), 'color': 'black'},
                           yshift=20)

    fig.update_xaxes(tickangle=-45)

    return fig


app.run_server(debug=True)
