# hsv_simulation

## 准备步骤

1. 安装`Microsoft Visual Studio`与`matlab`
2. 在`Microsoft Visual Studio`中安装并配置[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)库
3. 参考[教程](https://zhuanlan.zhihu.com/p/656804199)配置libtorch库与visual studio以实现python训练模型嵌入仿真（如果不需要libtorch可以注释相关内容编译速度更快）

## 操作步骤

1. 使用Microsoft Visual Studio打开`hsv_simulation_cpp.sln`文件
2. 编译并执行项目（需要配置好'eigen'与'libtorch'库）
3. 在matlab中运行画图文件`draw_cpp.m`
