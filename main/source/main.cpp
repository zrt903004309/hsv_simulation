#include<Global.h>
#include<Record.h>
#include <Controller_State.h>
#include <Vehicle_State.h>
#include <Guidance_State.h>

void main()
{
	//loadAeroCoefficients(0.05, 0.05);			// �����Ҫ��������ϵ����ͼ��ע�����������ļ�ֻ���д˺���
	ModelConfig mconfig;						// ����ģ�Ͳ�����ʼ��

	// ���������״̬���ʼ��
	VehicleState Vehicle_State(mconfig.h, mconfig.v, mconfig.alpha, mconfig.beta, mconfig.mu);					
	ControllerState Controller_State;			// ���������״̬���ʼ��
	VehiclePara Vehicle_Para;					// ����������������ʼ��	
	GuidanceState Guidance_State;				// ������̬�ǳ�ʼ��

	// �������ʼ��
	torch::jit::script::Module module;
	// һ��Ҫ������.pt�ļ�
	std::string model_name = "E:\\code\\codelearning\\hsv_thesis\\models\\model_with_weights.pt";
	// ����ʹ��cpu���������Ҫ��cuda����⿪ע�ͣ�����ȷ������������Ҫ������cuda
	////torch::Device device(torch::kCUDA);
	if (mconfig.network_en == 1) {
		// ���Լ���ģ��
		module = torch::jit::load(model_name);
		//module = torch::jit::load(model_name, device);
		std::cout << "ģ�ͼ��سɹ���" << std::endl;
	}
	
	// ��ӡģ�ͽṹ
	//std::cout << "\nģ�ͽṹ:\n" << module.dump_to_str(true, true, true) << std::endl;

	// Ӳ��ƽ̨��ʼ��
	CSerialPort mySerialPort;
	if (mconfig.hardware_en == 1) {
		if (!mySerialPort.InitPort(3, CBR_115200, 'E', 8, 1, EV_RXCHAR)) // �Ƿ�򿪴��ڣ�3�������������ӵ��Ե�com�ڣ��������豸�������鿴��Ȼ������������
		{
			std::cout << "initPort fail !" << std::endl;
		}
		else
		{
			std::cout << "initPort success !" << std::endl;
		}
		if (!mySerialPort.OpenListenThread())//�Ƿ�򿪼����̣߳������߳��������䷵��ֵ
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
			std::cout << "�����ɹ���" << std::endl;
		}
		Sleep(10000);
		mySerialPort.WriteData(0x04, 120.0f, 8);						// �����Ƿ���Գɹ�����
	}

	// delta[0]��Ӧ���� - delta_e delta[1]��Ӧ���� - delta_a delta[2]��Ӧ����� - delta_r

	//��ʱ
	clock_t start, stop;
	double duration;
	
	// ������ڼ�¼��txt
	FILE *File_Vehicle = NULL;
	FILE *File_Control = NULL;

	// ��Ĭ��·����(Ҳ���޸�·��)����д�뷽ʽ��txt���ڼ�¼���ݣ��������ڸ�txt����Զ�����
	File_Vehicle = fopen(mconfig.vehicle_filename.c_str(), "w");
    File_Control = fopen(mconfig.control_filename.c_str(), "w");

	start = clock();  //��ʼ��ʱ
	for(unsigned int i = 1; i <= mconfig.iters; i++)
	{
		if(i % mconfig.tctrl == 1)
		{  // ���������¡�����Ϊ������״̬��������״̬�����������������Ϊ������״̬
			Guidance_State.Guidance_State_Update(mconfig.step * i, mconfig);
			//Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig, module);
			//Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig);
			//Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig, mySerialPort);
			Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig, module, mySerialPort);
			Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // ��¼�������ݣ����txt
		}
		// ������״̬���¡�����Ϊ������״̬��������״̬�����������������沽�������Ϊ������״̬
		//Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // ��¼�������ݣ����txt
		Vehicle_State.Vehicle_State_Update(Controller_State, Vehicle_Para, mconfig.step, mconfig);
	}
	stop = clock();  // ������ʱ
	duration = (double)(stop - start) / CLK_TCK;  // ��Ϊ��ʱ���Ժ���Ϊ��λ�ģ�����CLK_TCK���Ի�����룬���Ϊ����ѭ�������ʱ��/s
	std::cout << duration << std::endl;
	printf("%f ��\n",duration);  // �������ѭ�������ʱ��

	if (mconfig.hardware_en == 1) {
		if (mySerialPort.WriteData(0x00, 1.0f, 8)) {
			std::cout << "�ػ��ɹ���" << std::endl;
		}
	}
	
	//// ����ģ������
	//torch::Tensor input = torch::tensor(
	//	{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 },
	//	torch::kFloat64
	//).view({ 1, 10 });  // ������״Ϊ[1, 10]��batch_size=1���ɼ�.to(device)
	////std::cout << "����������״: " << input.sizes() << std::endl;
	//std::vector<torch::jit::IValue> inputs;
	//inputs.push_back(input);


	//// ǰ������
	//at::Tensor output = module.forward(inputs).toTensor();
	////std::cout << "���������״: " << output.sizes() << std::endl;
	////std::cout << "���ֵ: " << output << std::endl;
	//// ��ȡ���
	//double* data = output.data_ptr<double>();
	//for (int64_t i = 0; i < 4; ++i) {
	//	std::cout << data[i] << " ";
	//}

}