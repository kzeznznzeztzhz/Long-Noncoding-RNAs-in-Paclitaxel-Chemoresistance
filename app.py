import pandas as pd
import numpy as np
import networkx as nx
import re
import plotly.graph_objects as go
import dash
from dash import html, dcc, Input, Output
import dash_bootstrap_components as dbc

# Load data
raw_df = np.array(pd.read_excel('function_human.xls', sheet_name='Sheet2'))[1:, :]
filt_ID_bio_process = raw_df[:, [1, 3]]
raw_int_df = np.array(pd.read_excel('Heatmap_interpretation.xlsx', sheet_name='Sheet1'))[:, [0, 7]]

def build_gene_ontology_graph(filt_ID_bio_process, raw_int_df):
    G = nx.Graph()
    center_node = "PAX"
    G.add_node(center_node, layer=0)
    gene_regulation = {gene: status.strip() for gene, status in raw_int_df}

    for gene, processes_str in filt_ID_bio_process:
        processes = [p.strip().lower() for p in re.split(r',|\band\b', processes_str) if p.strip()]
        for process in processes:
            G.add_node(process, layer=1)
            G.add_node(gene, layer=2, regulation=gene_regulation.get(gene, 'Unknown'))
            G.add_edge(process, gene)
            G.add_edge(center_node, process)
    return G

def get_node_color(layer, regulation):
    if layer == 0:
        return 'lightgreen'
    elif layer == 1:
        return 'lightblue'
    elif regulation == 'upregulation':
        return 'red'
    elif regulation == 'downregulation':
        return 'blue'
    else:
        return 'gray'

def generate_figure(G, focus_node=None):
    pos = nx.spring_layout(G, k=0.5, iterations=100, seed=42)
    nodes_to_show = G.nodes()
    edges_to_show = G.edges()

    if focus_node:
        neighbors = list(G.neighbors(focus_node))
        nodes_to_show = [focus_node] + neighbors
        edges_to_show = [(focus_node, n) for n in neighbors if G.has_edge(focus_node, n)]

    def curved_edges(u, v, curvature=0.2, resolution=10):
        # mid-point
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        # perpendicular offset
        dx, dy = x1 - x0, y1 - y0
        px, py = -dy, dx
        # normalize
        length = np.sqrt(px ** 2 + py ** 2)
        px, py = px / length, py / length
        cx, cy = mx + curvature * px, my + curvature * py
        # bezier curve points
        t = np.linspace(0, 1, resolution)
        x_curve = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * cx + t ** 2 * x1
        y_curve = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * cy + t ** 2 * y1
        return x_curve, y_curve

    edge_traces = []
    for u, v in edges_to_show:
        x_curve, y_curve = curved_edges(u, v, curvature=0.15, resolution=20)
        edge_traces.append(go.Scatter(
            x=x_curve, y=y_curve,
            mode='lines',
            line=dict(width=1.5, color='gray'),
            hoverinfo='none',
            showlegend=False
        ))

    node_x, node_y, node_text, node_color = [], [], [], []
    for node in nodes_to_show:
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(str(node))
        layer = G.nodes[node].get('layer')
        regulation = G.nodes[node].get('regulation', '').lower()
        node_color.append(get_node_color(layer, regulation))

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=node_text,
        textposition='top center',
        hoverinfo='text',
        marker=dict(size=20, color=node_color, line=dict(width=2, color='black')),
        name='Nodes',
        showlegend=False
    )

    # Legend
    legend_traces = [
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=15, color='lightgreen'),
                   name='PAX (center)'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=15, color='lightblue'),
                   name='Biological Process'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=15, color='red'),
                   name='Upregulated Gene'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=15, color='blue'),
                   name='Downregulated Gene'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=15, color='gray'),
                   name='Unknown Regulation'),
    ]

    fig = go.Figure(data=edge_traces + [node_trace] + legend_traces,
                    layout=go.Layout(
                        title='Gene Ontology Graph (Curved Edges + Click to Focus)',
                        showlegend=True,
                        legend=dict(itemsizing='constant', font=dict(size=12)),
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, visible=False),
                        yaxis=dict(showgrid=False, zeroline=False, visible=False),
                        height=800
                    ))

    return fig


# Build graph
G = build_gene_ontology_graph(filt_ID_bio_process, raw_int_df)

# Dash App
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "Gene Ontology Viewer"

app.layout = dbc.Container([
    html.H1("Gene Ontology Graph (Click a Node to Focus)", className='text-center my-3'),
    dcc.Graph(id='ontology-graph', figure=generate_figure(G), config={'scrollZoom': True}),
    dbc.Button("Reset Graph", id='reset-button', color='secondary', className='mt-3')
], fluid=True)

@app.callback(
    Output('ontology-graph', 'figure'),
    Input('ontology-graph', 'clickData'),
    Input('reset-button', 'n_clicks'),
    prevent_initial_call=True
)
def update_graph(click_data, reset_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return generate_figure(G)
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if trigger_id == 'reset-button':
        return generate_figure(G)
    elif click_data:
        clicked_node = click_data['points'][0]['text']
        return generate_figure(G, focus_node=clicked_node)

    return generate_figure(G)

# Run App
server = app.server
if __name__ == '__main__':

    app.run(debug=True)
