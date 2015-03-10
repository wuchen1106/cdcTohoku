#!/bin/bash

#nohup ./myminuit 189 1 def 0 -845 1.121 1.131 &
#nohup ./myminuit 193 1 def 0 -845 1.136 1.148 &
#nohup ./myminuit 194 1 def 0 -845 1.136 1.148 &
#nohup ./myminuit 218 1 def 0 -835 1.120 1.126 &
#nohup ./myminuit 221 1 def 0 -835 1.120 1.126 &
#nohup ./myminuit 222 1 def 0 -835 1.120 1.126 &
#nohup ./myminuit 250 1 flip 0 -835 1.082 1.097 &
#nohup ./myminuit 250 1 def2 0 -835 1.082 1.097 &
#nohup ./myminuit 250 1 def 0 -835 1.082 1.097 &
#nohup ./myminuit 251 1 def2 0 -838 1.05 1.05 -838 -836 -837 &
#nohup ./myminuit 251 1 fit 0 -838 1 1 -838 -836 -837 &
#nohup ./myminuit 251 1 flip 100000 -835 1.082 1.097 &

##./getFitGraphs 251 1 ori 0 0.229 2.470 0.313 4.116
#nohup ./getFitGraphs 251 1 def 0 0.229 2.470 0.313 4.116 &
#nohup ./getFitGraphs 251 1 def2 0 0.229 2.470 0.313 4.116 &
#nohup ./getFitGraphs 250 1 deflipf 0 0.165 1.255 0.178 1.883 &
#nohup ./getFitGraphs 250 1 flip 0 0.165 1.255 0.178 1.883 &
##./getFitGraphs 222 1 def 0 0.165 1.255 0.178 1.883
#nohup ./getFitGraphs 221 1 def 0 0.586 6.373 0.515 5.650 &

# ========================================

#nohup ./myminuit 250 1 it 0 0 500 &
#nohup ./myminuit 251 1 it 0 0 500 &

#nohup ./getFitGraphs 250 1 it 0 0 0.165 1.255 0.178 1.883 &
#nohup ./getFitGraphs 251 1 it 0 0 0.229 2.470 0.313 4.116 &

#nohup ./myminuit 250 1 it 0 1 200 &
#nohup ./myminuit 251 1 it 0 1 200 &

#nohup ./getFitGraphs 250 1 it 1 0 0.165 1.255 0.178 1.883 &
#nohup ./getFitGraphs 251 1 it 1 0 0.229 2.470 0.313 4.116 &

# ========================================

#nohup ./myminuit 142 1 nocut 0 0 100 &
#nohup ./myminuit 189 1 nocut 0 0 100 &
#nohup ./myminuit 193 1 nocut 0 0 100 &
#nohup ./myminuit 194 1 nocut 0 0 100 &
#nohup ./myminuit 215 1 nocut 0 0 100 &
#nohup ./myminuit 216 1 nocut 0 0 100 &
#nohup ./myminuit 217 1 nocut 0 0 100 &
#nohup ./myminuit 218 1 nocut 0 0 100 &
#nohup ./myminuit 219 1 nocut 0 0 100 &
#nohup ./myminuit 220 1 nocut 0 0 100 &
#nohup ./myminuit 221 1 nocut 0 0 100 &
#nohup ./myminuit 222 1 nocut 0 0 100 &
#nohup ./myminuit 250 1 nocut 0 0 100 &
#nohup ./myminuit 251 1 nocut 0 0 100 &

#nohup ./getFitGraphs 142 1 nocut 0 0 7 &
#nohup ./getFitGraphs 189 1 nocut 0 0 7 &
#nohup ./getFitGraphs 193 1 nocut 0 0 7 &
#nohup ./getFitGraphs 194 1 nocut 0 0 7 &
#nohup ./getFitGraphs 215 1 nocut 0 0 7 &
#nohup ./getFitGraphs 216 1 nocut 0 0 7 &
#nohup ./getFitGraphs 217 1 nocut 0 0 7 &
#nohup ./getFitGraphs 218 1 nocut 0 0 7 &
#nohup ./getFitGraphs 219 1 nocut 0 0 7 &
#nohup ./getFitGraphs 220 1 nocut 0 0 7 &
#nohup ./getFitGraphs 221 1 nocut 0 0 7 &
#nohup ./getFitGraphs 222 1 nocut 0 0 7 &
#nohup ./getFitGraphs 250 1 nocut 0 0 7 &
#nohup ./getFitGraphs 251 1 nocut 0 0 7 &

#nohup ./myminuit 142 1 nocut 0 1 100 &
#nohup ./myminuit 189 1 nocut 0 1 100 &
#nohup ./myminuit 221 1 nocut 0 1 100 &
#nohup ./myminuit 222 1 nocut 0 1 100 &
#nohup ./myminuit 250 1 nocut 0 1 100 &
#nohup ./myminuit 251 1 nocut 0 1 100 &

#nohup ./getFitGraphs 142 1 nocut 1 0 7 &
#nohup ./getFitGraphs 189 1 nocut 1 0 7 &
#nohup ./getFitGraphs 221 1 nocut 1 0 7 &
#nohup ./getFitGraphs 222 1 nocut 1 0 7 &
#nohup ./getFitGraphs 250 1 nocut 1 0 7 &
#nohup ./getFitGraphs 251 1 nocut 1 0 7 &

#nohup ./getFitGraphs 221 1 nocut 3 0 7 &

# ========================================

#nohup ./myminuit 216 1 hcut 0 0 100 &
#nohup ./myminuit 221 1 hcut 0 0 100 &
#nohup ./myminuit 222 1 hcut 0 0 100 &
#nohup ./myminuit 251 1 hcut 0 0 100 &

#nohup ./getFitGraphs 216 1 hcut 0 0 7 &
#nohup ./getFitGraphs 221 1 hcut 0 0 7 &
#nohup ./getFitGraphs 222 1 hcut 0 0 7 &
#nohup ./getFitGraphs 251 1 hcut 0 0 7 &

#nohup ./myminuit 221 1 hcut 0 1 100 &

#nohup ./getFitGraphs 221 1 hcut 1 0 7 &

#nohup ./myminuit 221 1 hcut 0 2 100 &

#nohup ./getFitGraphs 221 1 hcut 2 0 7 &

# ========================================
#nohup ./myminuit 221 1 t0 0 0 100 &
#nohup ./myminuit 251 1 t0 0 0 100 &

#nohup ./getFitGraphs 221 1 t0 0 0 7 &
#nohup ./getFitGraphs 251 1 t0 0 0 7 &

# ========================================
#nohup ./myminuit 142 1 0309 0 0 100 &
#nohup ./myminuit 189 1 0309 0 0 100 &
#nohup ./myminuit 193 1 0309 0 0 100 &
#nohup ./myminuit 194 1 0309 0 0 100 &
#nohup ./myminuit 215 1 0309 0 0 100 &
#nohup ./myminuit 216 1 0309 0 0 100 &
#nohup ./myminuit 217 1 0309 0 0 100 &
#nohup ./myminuit 218 1 0309 0 0 100 &
#nohup ./myminuit 219 1 0309 0 0 100 &
#nohup ./myminuit 220 1 0309 0 0 100 &
#nohup ./myminuit 221 1 0309 0 0 100 &
#nohup ./myminuit 222 1 0309 0 0 100 &
#nohup ./myminuit 250 1 0309 0 0 100 &
#nohup ./myminuit 250 1 hcut0 0 0 100 &
#nohup ./myminuit 250 1 hcut2 0 0 100 &
#nohup ./myminuit 251 1 0309 0 0 100 &

#nohup ./getFitGraphs 142 1 0309 0 0 7 &
nohup ./getFitGraphs 189 1 0309 0 0 7 &
nohup ./getFitGraphs 193 1 0309 0 0 7 &
nohup ./getFitGraphs 194 1 0309 0 0 7 &
#nohup ./getFitGraphs 215 1 0309 0 0 7 &
#nohup ./getFitGraphs 216 1 0309 0 0 7 &
#nohup ./getFitGraphs 217 1 0309 0 0 7 &
#nohup ./getFitGraphs 218 1 0309 0 0 7 &
#nohup ./getFitGraphs 219 1 0309 0 0 7 &
#nohup ./getFitGraphs 220 1 0309 0 0 7 &
#nohup ./getFitGraphs 221 1 0309 0 0 7 &
#nohup ./getFitGraphs 222 1 0309 0 0 7 &
#nohup ./getFitGraphs 250 1 0309 0 0 7 &
#nohup ./getFitGraphs 250 1 hcut0 0 0 7 &
#nohup ./getFitGraphs 250 1 hcut2 0 0 7 &
#nohup ./getFitGraphs 251 1 0309 0 0 7 &

#nohup ./myminuit 142 1 0309 0 1 100 &
#nohup ./myminuit 189 1 0309 0 1 100 &
#nohup ./myminuit 221 1 0309 0 1 100 &
#nohup ./myminuit 222 1 0309 0 1 100 &
#nohup ./myminuit 250 1 0309 0 1 100 &
#nohup ./myminuit 251 1 0309 0 1 100 &

#nohup ./getFitGraphs 142 1 0309 1 0 7 &
#nohup ./getFitGraphs 189 1 0309 1 0 7 &
#nohup ./getFitGraphs 221 1 0309 1 0 7 &
#nohup ./getFitGraphs 222 1 0309 1 0 7 &
#nohup ./getFitGraphs 250 1 0309 1 0 7 &
#nohup ./getFitGraphs 251 1 0309 1 0 7 &

#nohup ./getFitGraphs 221 1 0309 3 0 7 &
