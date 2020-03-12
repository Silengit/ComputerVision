#第二题运行指令
##简单版： 
```matlab
img1 = myImageFilter(img0);
%默认核为4×4的平均滤波器，默认步长水平与竖直都为1
```

##进阶版： 
```matlab
img1 = myImageFilter(img0, k);
%使用核k(二维矩阵)作为滤波器，默认步长水平与竖直都为1
```

##复杂版： 
```matlab
img1 = myImageFilter(img0, k, s);
%使用核k(二维矩阵)作为滤波器，使用s（1×2的向量，其中分量1为竖直步长，分量2为水平步长）作为步长
```

#第三题第一问运行指令
```matlab
img2 = globalHistogramEqualization(img0);
```

#第三题第一问运行指令
```matlab
img3 = localHistogramEqualization(img0);
```