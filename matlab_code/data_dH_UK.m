% hospital admissions data for UK available at 
% https://coronavirus.data.gov.uk/details/healthcare

data_dH_GB_inv2=[289  %  
338
343
363
359
356
361
379
447
443
440
498
426
461
481
544
568
579
576
537
637
597
704
713
832
822
770
740
862
970
1021
1120
1121
1154
1085
1268
1334
1407
1489
1501
1535
1426
1506
1676
1724
1753
1920
1849
1891
2003
2143
2355
2479
2600
2401
2471
2683
3041
3058
3144
2930
3110
3112
3391
3632
4044
4158
4063
3748
3938
3675
4099
4223
4577
4310
3931
4125
3997
4291
4394
4117
3937
3649
3439
3369
2916
3285
3250
3131
2871
2529
2183
2082
2274
2478
2368
2167
2044
1920
2144
2078
1957
1953
1753
1799
1686
1751
1753
1784
1658
1470
1378
1433
1502
1489
1490
1404
1384
1221
1280
1423
1613
1496
1587
1488
1460
1547
1685
1782
1715
1768
1670
1612
1647
1914
1971
1848
1763
1687
1555
1414
1586
1677
1514
1566
1488
1459
1373
1584
1485
1488
1532
1420
1256
1289
1275
1259
1186
1094
1086
1025
911
1044
980
1021
884
823
720
703
710
683
722
623
628
532
528
523
452
440
443
340
351
383
422
421
350
369
289
260
262
264
284
271
255
226
216
231
219
202
184
142
138
130
153
145
146
134
91
77
110
137
152
132
121
91
74
72
109
98
115
130
96
80
99
102
130
101
114
86
106
141
117
101
140
127
110
93
128
114
132
132
119
130
109
119
158
200
158
126
130
126
105
144
166
166
185
189
114
177
176
213
208
214
204
196
202
214
218
269
295
220
208
287
340
324
368
363
394
288
257
336
403
432
455
440
364
385
392
446
512
492
412
407
406
475
565
633
668
604
528
568
614
702
712
698
602
587
643
739
842
851
881
848
792
813
870
909
966
983
1029
915
943
910
1165
1206
1146
1250
1188
1180
1262
1443
1489
1517
1538
1355
1255
1525
1513
1527
1713
1737
1657
1629
1966
2012
1880
2088
2160
2066
2162
2463
2629
2967
3150
3046
2916
2921
2930
3284
3565
2832
3190
2489
2191
2227
1929
2085
1720
1273];

data_dH2=data_dH_GB_inv2(end:-1:1); % from 23 Mar

cum_sum_H0=4876-data_dH2(1);
data_cumsum_H2=[cum_sum_H0; cum_sum_H0+cumsum(data_dH2)];