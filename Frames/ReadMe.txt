{\rtf1\ansi\ansicpg1252\cocoartf1348\cocoasubrtf170
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural

\f0\fs24 \cf0 The file \'93Array_Parity_Check_Matrix.m\'94 includes the implementation of the binary measurement matrices constructed based on parity check matrices of Array codes. \
This construction receives two inputs , \'93n\'94 and \'93j\'94 which show the dimensionality of the unknown vector \'93x\'94 (or the number of the columns in the measurement matrix) and the column-weight respectively and outputs a binary matrix of \'93jq*q^2\'94 dimensions. \
This matrix is proved to have girth equal to 6. \
\
In our paper [1], we have proved that binary matrices with fixed column-weight are able to satisfy Robust Null-Space Property (RNSP) and achieve robust sparse recovery. We have also proved that in order to get the highest compression rate in any compressive sensing problem (using RNSP as the recovery guarantee), matrices with girth 6 are enough!\
\
-If you use this code in your research please support us by citing our paper in [1].\
\
[1] Compressed Sensing Using Binary Matrices of Nearly Optimal Dimensions, Mahsa Lotfi and Mathukumalli Vidyasagar, arXiv: 1808.03001\
\
   }