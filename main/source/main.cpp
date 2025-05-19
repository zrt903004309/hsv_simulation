#include<Global.h>
#include<Record.h>
#include <Controller_State.h>
#include <Vehicle_State.h>
#include <Guidance_State.h>

void main()
{
	//loadAeroCoefficients(0.05, 0.05);			// 如果需要产生气动系数画图，注释其他所有文件只运行此函数
	ModelConfig mconfig;						// 仿真模型参数初始化

	// 定义飞行器状态与初始化
	VehicleState Vehicle_State(mconfig.h, mconfig.v, mconfig.alpha, mconfig.beta, mconfig.mu);					
	ControllerState Controller_State;			// 定义控制器状态与初始化
	VehiclePara Vehicle_Para;					// 定义飞行器参数与初始化	
	GuidanceState Guidance_State;				// 跟踪姿态角初始化

	// 神经网络初始化
	torch::jit::script::Module module;
	// 一定要导出成.pt文件
	std::string model_name = "E:\\code\\codelearning\\hsv_thesis\\models\\model_with_weights.pt";
	// 现在使用cpu推理，如果需要用cuda推理解开注释，并且确保所有运算需要搬运至cuda
	////torch::Device device(torch::kCUDA);
	if (mconfig.network_en == 1) {
		// 尝试加载模型
		module = torch::jit::load(model_name);
		//module = torch::jit::load(model_name, device);
		std::cout << "模型加载成功！" << std::endl;
	}
	
	// 打印模型结构
	//std::cout << "\n模型结构:\n" << module.dump_to_str(true, true, true) << std::endl;

	// 硬件平台初始化
	CSerialPort mySerialPort;
	if (mconfig.hardware_en == 1) {
		if (!mySerialPort.InitPort(3, CBR_115200, 'E', 8, 1, EV_RXCHAR)) // 是否打开串口，3就是你外设连接电脑的com口，可以在设备管理器查看，然后更改这个参数
		{
			std::cout << "initPort fail !" << std::endl;
		}
		else
		{
			std::cout << "initPort success !" << std::endl;
		}
		if (!mySerialPort.OpenListenThread())//是否打开监听线程，开启线程用来传输返回值
		{
			std::cout << "OpenListenThread fail !" << std::endl;
		}
		else
		{
			std::cout << "OpenListenThread success !" << std::endl;
		}
		float data = 0.0f;
		unsigned char control = 0x01;
		if (mySerialPort.WriteData(control, data, 8)) {
			std::cout << "开机成功！" << std::endl;
		}
		Sleep(10000);
		mySerialPort.WriteData(0x04, 120.0f, 8);						// 测试是否可以成功运行
	}

	// delta[0]对应左翼 - delta_e delta[1]对应右翼 - delta_a delta[2]对应方向舵 - delta_r

	//计时
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
			//Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig, module);
			//Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig);
			//Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig, mySerialPort);
			Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig, module, mySerialPort);
			Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // 记录所得数据，输出txt
		}
		// 飞行器状态更新。输入为飞行器状态、控制器状态、飞行器参数、仿真步长。输出为飞行器状态
		//Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // 记录所得数据，输出txt
		Vehicle_State.Vehicle_State_Update(Controller_State, Vehicle_Para, mconfig.step, mconfig);
	}
	stop = clock();  // 结束计时
	duration = (double)(stop - start) / CLK_TCK;  // 因为计时是以毫秒为单位的，除以CLK_TCK可以换算成秒，输出为迭代循环所需的时间/s
	std::cout << duration << std::endl;
	printf("%f 秒\n",duration);  // 输出迭代循环所需的时间

	if (mconfig.hardware_en == 1) {
		if (mySerialPort.WriteData(0x00, 1.0f, 8)) {
			std::cout << "关机成功！" << std::endl;
		}
	}
	
	//// 创建模型输入
	//torch::Tensor input = torch::tensor(
	//	{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 },
	//	torch::kFloat64
	//).view({ 1, 10 });  // 调整形状为[1, 10]，batch_size=1，可加.to(device)
	////std::cout << "输入张量形状: " << input.sizes() << std::endl;
	//std::vector<torch::jit::IValue> inputs;
	//inputs.push_back(input);


	//// 前向推理
	//at::Tensor output = module.forward(inputs).toTensor();
	////std::cout << "输出张量形状: " << output.sizes() << std::endl;
	////std::cout << "输出值: " << output << std::endl;
	//// 读取输出
	//double* data = output.data_ptr<double>();
	//for (int64_t i = 0; i < 4; ++i) {
	//	std::cout << data[i] << " ";
	//}

}