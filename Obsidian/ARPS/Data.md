| Datat type | time | time interval | index |
| --- | --- | --- | --- |
| ar2016062307 | 0700-0728 UTC | 1 second | 0-435, 436-822, 823-899, 900-1380, 1381-1680, 1202 
| ar2016062306 | 0630-0700 UTC | 5 seconds | 0-1795 | 
| ar2016062305 | 0530-0630 UTC | 5 seconds | 0-3600 |


| Data type | Data dir | time interval |
| --- | --- | --- |
| ar2016062307 | /data/data_s12_1/whuang/arpsdiag/ | 1 second |
| ar2016062307 | /gpfs3/whuang/arps/cut_for_arps/arpsdata_0000_2857/ | 3 seconds |
| ar2016062306 | /gpfs3/whuang/arps/cut_for_arps/arpsdata_201606230630/ | 5 seconds |
| ar2016062305 | /gpfs3/whuang/arps/cut_for_arps/arpsdata_201606230530/ | 5 seconds |


| Data type | solv_opt | Diag dir | time range |
| --- | --- | --- | --- |
| ar2016062307 | 4 | /data/data_s12_1/whuang/arpsdiag/ | 0-1680 (1 second) |
| ar2016062307 | 4 | /gpfs3/whuang/arps/cut_for_arps/arpsdata_0000_2857/ <br> *__the total dir of ==ar206062307==, the data in /data/data_s12_1/whuang/arpsdiag/ has linked here__* |
| ar2016062307 | 6 | /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0700/ | 0-1680 (5 and 10 seconds) |
| ar2016062306 | 6 | /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0630/ | 0-1790 (10 seconds) |
| ar2016062305 | 6 | /data/data_s13_3/whuang/arpsdiag/ | 600-3395 (5 seconds) | 
| ar2016062305 | 6 | /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0530/ <br> *__the total dir of include ==all files for solv_opt == 6== (ar2016062305, ar2016062306, and ar2016062307), all files have linked here__* | 0-7800  | 

PS. The dir /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0530/ include many files in different directories. 
- /data/data_s13_3/whuang/arpsdiag/
- /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0729/
- /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0700/
- /gpfs3/whuang/arps/arpsdiag/forHuangWei/diag_output/test/0630/