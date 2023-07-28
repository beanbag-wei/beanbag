### ARPSTRAJC
1. 改动了arpstrajc和arpscalcgtrajc模块，此前是依靠ARPS里的时间数据辨别时间序列，如果ARPS数据是断开的那么就会出现问题（直接run arps不会出现这个问题，但是wrf2arps中途很有可能断掉，那么arps内部时间就会出现问题），现在的版本是按hdf文件的时间序列读取的，可以规避这个问题
2.  ARPSDIAG其实也有上述问题，目前先不管 
3. 之前计算的版本/xue_home/whuang/gpfs3/arps5.4.1编译存在问题，目前使用/xue_home/whuang/gpfs3/arps5.4.1_ynan版本
***
### ARPSDIAG
1. solv_opt == 4计算得到的扰动气压是有问题的。因为solv_opt == 4计算前重新计算了base state的气压，而这个计算得到的扰动气压是不对的，气压梯度力不连续，目前采取不重新计算base state的气压值，还是使用hdfgrdbas的平均气压场作为平均气压，从而得到扰动气压
2. 
