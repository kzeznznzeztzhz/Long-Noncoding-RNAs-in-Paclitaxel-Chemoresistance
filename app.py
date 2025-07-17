import dash
from dash import html, dcc, Input, Output
import dash_bootstrap_components as dbc
import os
import pandas as pd
import numpy as np
import networkx as nx
import re
import plotly.graph_objects as go
from flask import send_from_directory

# Load your data and build the graph
df1 = pd.read_excel('function_human.xls', sheet_name='Sheet2')
df2 = pd.read_excel('Heatmap_interpretation.xlsx', sheet_name='Sheet1')
raw_df = np.array(df1)[1:, :]
filt_ID_bio_process = raw_df[:, [1, 3]]
raw_int_df = np.array(df2)[:, [0, 7]]

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

    def curved_edges(u, v, curvature=0.15, resolution=20):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = x1 - x0, y1 - y0
        px, py = -dy, dx
        length = np.sqrt(px ** 2 + py ** 2)
        px, py = px / length, py / length
        cx, cy = mx + curvature * px, my + curvature * py
        t = np.linspace(0, 1, resolution)
        x_curve = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * cx + t ** 2 * x1
        y_curve = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * cy + t ** 2 * y1
        return x_curve, y_curve

    edge_traces = [
        go.Scatter(
            x=curved_edges(u, v)[0],
            y=curved_edges(u, v)[1],
            mode='lines',
            line=dict(width=1.5, color='gray'),
            hoverinfo='none',
            showlegend=False
        ) for u, v in edges_to_show
    ]

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
        showlegend=False
    )

    legend_traces = [
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=15, color=c), name=n)
        for c, n in [
            ('lightgreen', 'PAX (center)'),
            ('lightblue', 'Biological Process'),
            ('red', 'Upregulated Gene'),
            ('blue', 'Downregulated Gene'),
            ('gray', 'Unknown Regulation')
        ]
    ]

    fig = go.Figure(data=edge_traces + [node_trace] + legend_traces)
    fig.update_layout(
        title='Gene Ontology Graph (Tap a Node to Focus)',
        showlegend=True,
        hovermode='closest',
        margin=dict(b=80, l=10, r=10, t=40),
        height=800,
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.25,
            xanchor="center",
            x=0.5,
            font=dict(size=12)
        )
    )
    return fig

G = build_gene_ontology_graph(filt_ID_bio_process, raw_int_df)

app = dash.Dash(__name__, use_pages=False, suppress_callback_exceptions=True, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

app.layout = html.Div([
    dcc.Location(id='url'),
    html.Div(id='page-content')
])

@app.callback(Output('page-content', 'children'), Input('url', 'pathname'))
def display_page(pathname):
    if pathname == "/":
        return html.Div([
            html.H2("Appendix of Long Noncoding RNAs in Paclitaxel Chemoresistance", className='text-center my-4'),
            dbc.Button("Interactive Gene Ontology Graph", href="/app", color="primary", className="me-2"),
            dbc.Button("Detailed LINC-Cancer Documentation", href="/html", color="secondary")
        ], className='text-center')
    elif pathname == "/app":
        return dbc.Container([
            html.H1("Gene Ontology Graph (Tap a Node to Focus)", className='text-center my-3'),
            dbc.Row([
                dcc.Graph(
                    id='ontology-graph',
                    figure=generate_figure(G),
                    config={'scrollZoom': True},
                    style={"width": "100%", "height": "auto"}
                )
            ], className="gx-0"),
            dbc.Button("Reset Graph", id='reset-button', color='secondary', className='my-3')
        ], fluid=True, className="px-2")
    elif pathname == "/html":
        with open("cancer_st.html", 'r', encoding='utf-8') as f:
            html_content = f.read()
        return html.Iframe(srcDoc=html_content, style={"width": "100%", "height": "1000px", "border": "none"})
    return html.Div("404 - Page not found")

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

@app.server.route('/Poster_References.pdf')
def serve_pdf():
    return send_from_directory(os.path.join(app.root_path, 'static'), 'Poster_References.pdf')

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 8050))
    app.run(host="0.0.0.0", port=port)
