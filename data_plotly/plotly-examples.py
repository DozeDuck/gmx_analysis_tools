#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 20:40:15 2023

@author: dozeduck
"""

# xaxis
import plotly.graph_objs as go

x = [-2, -1, 0, 1, 2]
y = [-1, 1, 0, 2, -2]

trace = go.Scatter(x=x, y=y, mode='markers')
data = [trace]
layout = go.Layout(title='My Scatter Plot with Origin Superimposed',
                   xaxis=dict(range=[-2, 2], zeroline=True),
                   # xaxis=dict(autorange=True, zeroline=True),
                   yaxis=dict(range=[-2, 2], zeroline=True))

fig = go.Figure(data=data, layout=layout)
fig.show()


# Scatter Plot with Labels:
import plotly.graph_objs as go

x = [1, 2, 3, 4]
y = [10, 20, 30, 40]
labels = ['A', 'B', 'C', 'D']

trace = go.Scatter(x=x, y=y, mode='markers+text', text=labels, textposition='bottom center')
data = [trace]
layout = go.Layout(title='My Scatter Plot with Labels')
fig = go.Figure(data=data, layout=layout)
fig.show()

# Bar Chart:
import plotly.graph_objs as go

x = ['Apples', 'Bananas', 'Grapes']
y = [10, 5, 20]

trace = go.Bar(x=x, y=y, marker=dict(color='rgb(158,202,225)', line=dict(color='rgb(8,48,107)', width=1.5)))
data = [trace]
layout = go.Layout(title='My Bar Chart')
fig = go.Figure(data=data, layout=layout)
fig.show()

# Line Chart with Multiple Lines:
import plotly.graph_objs as go

x = [1, 2, 3, 4]
y1 = [10, 20, 15, 25]
y2 = [5, 15, 10, 20]

trace1 = go.Scatter(x=x, y=y1, name='Line 1')
trace2 = go.Scatter(x=x, y=y2, name='Line 2')

data = [trace1, trace2]
layout = go.Layout(title='My Line Chart with Multiple Lines')
fig = go.Figure(data=data, layout=layout)
fig.show()

# Pie Chart:
import plotly.graph_objs as go

labels = ['Apples', 'Bananas', 'Grapes']
values = [10, 5, 20]

trace = go.Pie(labels=labels, values=values)
data = [trace]
layout = go.Layout(title='My Pie Chart')
fig = go.Figure(data=data, layout=layout)
fig.show()


# Heatmap:
import plotly.graph_objs as go
import numpy as np

z = np.random.rand(10, 10)

trace = go.Heatmap(z=z)
data = [trace]
layout = go.Layout(title='My Heatmap')
fig = go.Figure(data=data, layout=layout)
fig.show()

