#include<Global.h>
#include<Record.h>
#include <Controller_State.h>
#include <Vehicle_State.h>
#include <Guidance_State.h>

using namespace std;

void main()
{
	//loadAeroCoefficients(0.05, 0.05);			// �����Ҫ��������ϵ����ͼ��ע�����������ļ�ֻ���д˺���
	ModelConfig mconfig;						// ����ģ�Ͳ�����ʼ��

	// ���������״̬���ʼ��
	VehicleState Vehicle_State(mconfig.h, mconfig.v, mconfig.alpha, mconfig.beta, mconfig.mu);					
	ControllerState Controller_State;			// ���������״̬���ʼ��
	VehiclePara Vehicle_Para;					// ����������������ʼ��	
	GuidanceState Guidance_State;				// ������̬�ǳ�ʼ��

	// delta[0]��Ӧ���� delta[1]��Ӧ���� delta[2]��Ӧ�����

	// ��ʱ
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
			Controller_State.Controller_State_Update(Vehicle_State, Guidance_State, Vehicle_Para, mconfig);
			Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // ��¼�������ݣ����txt
		}
		// ������״̬���¡�����Ϊ������״̬��������״̬�����������������沽�������Ϊ������״̬
		//Record(Vehicle_State, Controller_State, Guidance_State, File_Vehicle, File_Control);  // ��¼�������ݣ����txt
		Vehicle_State.Vehicle_State_Update(Controller_State, Vehicle_Para, mconfig.step);
	}
	stop = clock();  // ������ʱ
	duration = (double)(stop - start) / CLK_TCK;  // ��Ϊ��ʱ���Ժ���Ϊ��λ�ģ�����CLK_TCK���Ի�����룬���Ϊ����ѭ�������ʱ��/s
	cout << duration << endl;
	printf("%f ��\n",duration);  // �������ѭ�������ʱ��

}