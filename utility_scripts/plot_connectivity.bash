#!/bin/bash


#1 parser.add_argument("trackingdir", type=str)
#2 parser.add_argument("annualonlyswitch", type=bool)
#3 parser.add_argument("ignorestagnantswitch", type=bool)
#4 parser.add_argument("logscaleswitch", type=bool)
#5 parser.add_argument("plotlineswitch", type=bool)


python /home/blaughli/tracking_project_v2/processing/plotting/pdfs/connectivity/plot_connectivity_allPLDs_withPaulMagic.py $1 $2 $3 $4 $5
