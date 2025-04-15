#include<Global.h>
#include<Record.h>
#include <Controller_State.h>
#include <Vehicle_State.h>
#include <Guidance_State.h>

using namespace std;

void main()
{
	//loadAeroCoefficients(0.05, 0.05);			// 如果需要产生气动系数画图，注释其他所有文件只运行此函数
	ModelConfig mconfig;						// 仿真模型参数初始化

	// 定义飞行器状态与初始化
	VehicleState Vehicle_State(mconfig.h, mconfig.v, mconfig.alpha, mconfig.beta, mconfig.mu);					
	ControllerState Controller_State;			// 定义控制器状态与初始化
	VehiclePara Vehicle_Para;					// 定义飞行器参数与初始化	
	GuidanceState Guidance_State;				// 跟踪姿态角初始化

	// delta[0]对应左翼 delta[1]对应右翼 delta[2]对应方向舵

	// 计时
	clock_t start, stop;
	double duration;
	
	// 清空用于记录的txt
	FILE *File_Vehicle = NULL;
	FILE *File_Control = NULL;

	// 在默认路径下(也可修改路径)，以写入方式打开txt用于记录数据，若不存在该txt则会自动生成
	File_Vehicle = fopen(mconfig.vehicle_filename.c_str(), "w");
    File_Control = fopen(mconfig.control_filename.c_str(), "w");

	start = clock();  //开始计时
	for(unsigned int i = 1; i <= mconfig.iters; i++)
	{
		if(i % mconfig.tctrl == 1)
		{  // 控制器更新。输入为飞行器状态、控制器状态、飞行器参数。输出为控制器状态
			Guidance_State.Guidance_State_Update(mconfig.step * i, mconfig);
			Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig);
			Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // 记录所得数据，输出txt
		}
		// 飞行器状态更新。输入为飞行器状态、控制器状态、飞行器参数、仿真步长。输出为飞行器状态
		//Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // 记录所得数据，输出txt
		Vehicle_State.Vehicle_State_Update(Controller_State, Vehicle_Para, mconfig.step);
	}
	stop = clock();  // 结束计时
	duration = (double)(stop - start) / CLK_TCK;  // 因为计时是以毫秒为单位的，除以CLK_TCK可以换算成秒，输出为迭代循环所需的时间/s
	cout << duration << endl;
	printf("%f 秒\n",duration);  // 输出迭代循环所需的时间

}