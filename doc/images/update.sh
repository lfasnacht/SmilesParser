#!/bin/sh
dot -Teps Vanillin-flat.dot |epstopdf --filter > Vanillin-flat.pdf
dot -Teps Vanillin-tree.dot |epstopdf --filter > Vanillin-tree.pdf
dot -Teps callgraph.dot |epstopdf --filter > callgraph.pdf
