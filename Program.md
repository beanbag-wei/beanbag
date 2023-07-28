### read the arpsdiag file
  The file has negative values, missing values = -99999
  Therefore, 用空格作为分隔符是行不通的，在出现missing value的时候会出现error
  ==solution： 用’\s+‘作为分隔符，去匹配空格==
  '\\s+'和'\\s*'的含义：
  ![[1690291206078.png]]

### print with no limit
  np.set_printoptions(threshold=np.inf)

