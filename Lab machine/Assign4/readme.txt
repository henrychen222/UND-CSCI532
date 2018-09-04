1.编译：nvcc kernel.cu -o ass4
2.执行nvprof ./ass4
3.查看nvprof打印出来的执行时间，具体解释google
4.改代码：srcHeight是输入矩阵高，srcWidth是输入矩阵宽，maskHeight是卷积核高度，maskWidth是卷积核宽度，按照作业要求设置，编译后查看nvprof打印的执行时间就行了
