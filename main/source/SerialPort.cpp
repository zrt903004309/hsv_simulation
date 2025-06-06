//#include "StdAfx.h"
#include "SerialPort.h"
#include <process.h>
#include <iostream>

/** 线程退出标志 */ 
bool CSerialPort::s_bExit = false;
/** 当串口无数据时,sleep至下次查询间隔的时间,单位:秒 */ 
const UINT SLEEP_TIME_INTERVAL = 5;

DWORD findFrameStart(unsigned char* buffer, DWORD bufferSize);

CSerialPort::CSerialPort(void)
: m_hListenThread(INVALID_HANDLE_VALUE)
{
	m_hComm = INVALID_HANDLE_VALUE;
	m_hListenThread = INVALID_HANDLE_VALUE;

	InitializeCriticalSection(&m_csCommunicationSync);
	delta = 180.0f;
}

CSerialPort::~CSerialPort(void)
{
	CloseListenTread();
	ClosePort();
	DeleteCriticalSection(&m_csCommunicationSync);
}

bool CSerialPort::InitPort( UINT portNo /*= 1*/,UINT baud /*= CBR_9600*/,char parity /*= 'N'*/,
						    UINT databits /*= 8*/, UINT stopsbits /*= 1*/,DWORD dwCommEvents /*= EV_RXCHAR*/ )
{

	/** 临时变量,将制定参数转化为字符串形式,以构造DCB结构 */ 
	char szDCBparam[50];
	sprintf_s(szDCBparam, "baud=%d parity=%c data=%d stop=%d", baud, parity, databits, stopsbits);

	/** 打开指定串口,该函数内部已经有临界区保护,上面请不要加保护 */ 
	if (!openPort(portNo))
	{
		return false;
	}

	/** 进入临界段 */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** 是否有错误发生 */ 
	BOOL bIsSuccess = TRUE;

    /** 在此可以设置输入输出的缓冲区大小,如果不设置,则系统会设置默认值.
	 *  自己设置缓冲区大小时,要注意设置稍大一些,避免缓冲区溢出
	 */
	if (bIsSuccess )
	{
		bIsSuccess = SetupComm(m_hComm,32768,32768);
	}

	/** 设置串口的超时时间,均设为0,表示不使用超时限制 */
	COMMTIMEOUTS  CommTimeouts;
	CommTimeouts.ReadIntervalTimeout         = 0;
	CommTimeouts.ReadTotalTimeoutMultiplier  = 0;
	CommTimeouts.ReadTotalTimeoutConstant    = 0;
	CommTimeouts.WriteTotalTimeoutMultiplier = 0;
	CommTimeouts.WriteTotalTimeoutConstant   = 0; 
	if ( bIsSuccess)
	{
		bIsSuccess = SetCommTimeouts(m_hComm, &CommTimeouts);
	}

	DCB  dcb;
	if ( bIsSuccess )
	{
		// 将ANSI字符串转换为UNICODE字符串
		DWORD dwNum = MultiByteToWideChar (CP_ACP, 0, szDCBparam, -1, NULL, 0);
		wchar_t *pwText = new wchar_t[dwNum] ;
		if (!MultiByteToWideChar (CP_ACP, 0, szDCBparam, -1, pwText, dwNum))
		{
			bIsSuccess = TRUE;
		}

		/** 获取当前串口配置参数,并且构造串口DCB参数 */ 
		bIsSuccess = GetCommState(m_hComm, &dcb) && BuildCommDCB(pwText, &dcb) ;
		/** 开启RTS flow控制 */ 
		dcb.fRtsControl = RTS_CONTROL_ENABLE; 

		/** 释放内存空间 */ 
		delete [] pwText;
	}

	if ( bIsSuccess )
	{
		/** 使用DCB参数配置串口状态 */ 
		bIsSuccess = SetCommState(m_hComm, &dcb);
	}
		
	/**  清空串口缓冲区 */
	PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_TXCLEAR | PURGE_RXABORT | PURGE_TXABORT);

	/** 离开临界段 */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return bIsSuccess==TRUE;
}

bool CSerialPort::InitPort( UINT portNo ,const LPDCB& plDCB )
{
	/** 打开指定串口,该函数内部已经有临界区保护,上面请不要加保护 */ 
	if (!openPort(portNo))
	{
		return false;
	}
	
	/** 进入临界段 */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** 配置串口参数 */ 
	if (!SetCommState(m_hComm, plDCB))
	{
		return false;
	}

	/**  清空串口缓冲区 */
	PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_TXCLEAR | PURGE_RXABORT | PURGE_TXABORT);

	/** 离开临界段 */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return true;
}

void CSerialPort::ClosePort()
{
	/** 如果有串口被打开，关闭它 */
	if( m_hComm != INVALID_HANDLE_VALUE )
	{
		CloseHandle( m_hComm );
		m_hComm = INVALID_HANDLE_VALUE;
	}
}

bool CSerialPort::openPort( UINT portNo )
{
	/** 进入临界段 */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** 把串口的编号转换为设备名 */ 
    char szPort[50];
	sprintf_s(szPort, "COM%d", portNo);

	/** 打开指定的串口 */ 
	m_hComm = CreateFileA(szPort,		                /** 设备名,COM1,COM2等 */ 
						 GENERIC_READ | GENERIC_WRITE,  /** 访问模式,可同时读写 */   
						 0,                             /** 共享模式,0表示不共享 */ 
					     NULL,							/** 安全性设置,一般使用NULL */ 
					     OPEN_EXISTING,					/** 该参数表示设备必须存在,否则创建失败 */ 
						 0,    
						 0);    

	/** 如果打开失败，释放资源并返回 */ 
	if (m_hComm == INVALID_HANDLE_VALUE)
	{
		std::cout << "串口打开失败" << std::endl;
		LeaveCriticalSection(&m_csCommunicationSync);
		return false;
	}
	std::cout << "串口打开成功" << std::endl;

	/** 退出临界区 */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return true;
}

bool CSerialPort::OpenListenThread()
{
	/** 检测线程是否已经开启了 */ 
	if (m_hListenThread != INVALID_HANDLE_VALUE)
	{
		/** 线程已经开启 */ 
		return false;
	}

	s_bExit = false;
	/** 线程ID */ 
	UINT threadId;
	/** 开启串口数据监听线程 */ 
	m_hListenThread = (HANDLE)_beginthreadex(NULL, 0, ListenThread, this, 0, &threadId);
	if (!m_hListenThread)
	{
		return false;
	}
	/** 设置线程的优先级,高于普通线程 */ 
	if (!SetThreadPriority(m_hListenThread, THREAD_PRIORITY_ABOVE_NORMAL))
	{
		return false;
	}

	return true;
}

bool CSerialPort::CloseListenTread()
{	
	if (m_hListenThread != INVALID_HANDLE_VALUE)
	{
		/** 通知线程退出 */ 
		s_bExit = true;

		/** 等待线程退出 */ 
		Sleep(10);

		/** 置线程句柄无效 */ 
		CloseHandle( m_hListenThread );
		m_hListenThread = INVALID_HANDLE_VALUE;
	}
	return true;
}

UINT CSerialPort::GetBytesInCOM()
{
	DWORD dwError = 0;	/** 错误码 */ 
	COMSTAT  comstat;   /** COMSTAT结构体,记录通信设备的状态信息 */ 
	memset(&comstat, 0, sizeof(COMSTAT));

	UINT BytesInQue = 0;
	/** 在调用ReadFile和WriteFile之前,通过本函数清除以前遗留的错误标志 */ 
	if ( ClearCommError(m_hComm, &dwError, &comstat) )
	{
		BytesInQue = comstat.cbInQue; /** 获取在输入缓冲区中的字节数 */ 
	}

	return BytesInQue;
}

UINT WINAPI CSerialPort::ListenThread( void* pParam )
{
	/** 得到本类的指针 */ 
	CSerialPort *pSerialPort = reinterpret_cast<CSerialPort*>(pParam);

	// 线程循环,轮询方式读取串口数据
	while (!pSerialPort->s_bExit) 
	{
		UINT BytesInQue = pSerialPort->GetBytesInCOM();
		/** 如果串口输入缓冲区中无数据,则休息一会再查询 */ 
		if ( BytesInQue == 0 )
		{
			Sleep(SLEEP_TIME_INTERVAL);
			continue;
		}

		/** 读取输入缓冲区中的数据并输出显示 */
		char cRecved = 0x00;
		do
		{
			cRecved = 0x00;
			if(pSerialPort->ReadChar(cRecved) == true)
			{
				//std::cout << cRecved ; 
				continue;
			}
		}while(--BytesInQue);
	}

	return 0;
}

bool CSerialPort::ReadChar( char &cRecved )
{
	BOOL  bResult     = TRUE;
	DWORD BytesRead   = 0;
	if(m_hComm == INVALID_HANDLE_VALUE)
	{
		return false;
	}

	/** 临界区保护 */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** 从缓冲区读取一个字节的数据 */ 
	bResult = ReadFile(m_hComm, &cRecved, 1, &BytesRead, NULL);
	if ((!bResult))
	{ 
		/** 获取错误码,可以根据该错误码查出错误原因 */ 
		DWORD dwError = GetLastError();

		/** 清空串口缓冲区 */ 
		PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT);
		LeaveCriticalSection(&m_csCommunicationSync);

		return false;
	}

	/** 离开临界区 */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return (BytesRead == 1);

}

bool CSerialPort::ReadByte(unsigned char* buffer, DWORD bytesToRead) {
	BOOL bResult = TRUE;
	DWORD bytesRead = 0;

	if (m_hComm == INVALID_HANDLE_VALUE) {
		return false;
	}

	// 临界区保护
	EnterCriticalSection(&m_csCommunicationSync);

	// 从缓冲区读取指定长度字节的数据
	bResult = ReadFile(m_hComm, buffer, bytesToRead, &bytesRead, NULL);

	/*for (DWORD i = 0; i < bytesRead; ++i) {
		printf("0x%02X ", buffer[i]);
	}
	printf("\n");*/

	if (!bResult) {
		// 获取错误码，可以根据该错误码查出错误原因
		DWORD dwError = GetLastError();

		std::cout << "Failed to read data, error code: " << dwError << std::endl;
		// 清空串口缓冲区
		PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT);
		LeaveCriticalSection(&m_csCommunicationSync);

		return false;
	}

	// 离开临界区
	LeaveCriticalSection(&m_csCommunicationSync);

	return (bytesRead == bytesToRead);
}
float convertToSingle(const unsigned char data[]) {
	// 提取 data 数组中的第 4 到第 7 个元素
	unsigned char byteArray[4] = { data[3], data[4], data[5], data[6] };

	// 将 byteArray 转换为 float 类型
	float* floatPtr = reinterpret_cast<float*>(byteArray);
	return *floatPtr;
}

bool CSerialPort::WriteData(unsigned char control, float data, unsigned int length )
{
	BOOL   bResult     = TRUE;
	DWORD  BytesToSend = 0;
	if(m_hComm == INVALID_HANDLE_VALUE)
	{
		std::cout << "Serial port is not opened" << std::endl;
		return false;
	}

	/** 临界区保护 */ 
	//EnterCriticalSection(&m_csCommunicationSync);

	// 控制位，告知更改dsp哪个变量
	unsigned char stateIndex = static_cast<unsigned char>(control);

	unsigned char* pData = new unsigned char[8];
	pData[0] = 0xAA;
	pData[1] = 0x55;
	pData[2] = stateIndex;
	// 将float类型的数据转换为字节数组
	unsigned char* floatBytes = reinterpret_cast<unsigned char*>(&data);
	for (int i = 0; i < sizeof(float); ++i) {
		pData[3 + i] = floatBytes[i];
	}
	pData[7] = 0x0D;

	// 打印发送的指令
	float a = convertToSingle(pData);
	std::cout << (a-180) << std::endl;

	/** 向缓冲区写入指定量的数据 */ 
	bResult = WriteFile(m_hComm, pData, length, &BytesToSend, NULL);
	if (!bResult)  
	{
		DWORD dwError = GetLastError();
		std::cout << "Failed to write data, error code: " << dwError << std::endl;
		/** 清空串口缓冲区 */ 
		PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT);
		//LeaveCriticalSection(&m_csCommunicationSync);

		return false;
	}

	/** 离开临界区 */ 
	//LeaveCriticalSection(&m_csCommunicationSync);

	return true;
}


double CSerialPort::GetDelta(double cmdDelta, double stime) {
	cmdDelta += 180.0;
	float fcmdDelta = (float)cmdDelta;
	WriteData(0x04, fcmdDelta, 8);
	Sleep(stime * 1000);
	unsigned char* Data = new unsigned char[8];

	// 清空串口接收缓冲区
	if (!PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT)) {
		LeaveCriticalSection(&m_csCommunicationSync);
		/*return false;*/
	}

	float res = readFloatData();
	//std::cout << res << std::endl;
	return (double)(res-180);
}

bool CSerialPort::findValidFrame(unsigned char* buffer, DWORD bufferSize, float& result) {
	// 至少需要8个字节才可能包含完整帧
	if (bufferSize < 8) return false;

	// 滑动窗口搜索帧头
	for (DWORD i = 0; i <= bufferSize - 8; i++) {
		// 检查是否找到帧头
		if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2] == 0x11) {
		//if ((buffer[i]&0xFF) == 0xAA && (buffer[i + 1]&0xFF) == 0x55 && (buffer[i + 2]&0xFF) == 0x11) {
			//if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2] == 0x11) {
			// 检查帧尾
			//if ((buffer[i + 7] & 0xFF) == 0x0D) {
			if (buffer[i + 7] == 0x0D) {
				for (DWORD j = i; j < 8; ++j) {
					printf("0x%02X ", buffer[j]);
				}
				printf("\n");
				// 提取浮点数数据
				result = convertToSingle(&buffer[i + 3]);
				// 返回true表示找到有效帧，并通过引用返回结果
				return true;
			}
		}
	}
	return false;
}

float CSerialPort::readFloatData() {
	const int BUFFER_SIZE = 256;
	unsigned char buffer[BUFFER_SIZE] = { 0 };
	DWORD bufferPos = 0;
	float result = 0.0f;

	// 清空串口接收缓冲区
	if (!PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT)) {
		LeaveCriticalSection(&m_csCommunicationSync);
		/*return false;*/
	}

	while (true) {
		 // 尝试读取更多数据
		if (bufferPos < BUFFER_SIZE - 8) {
			unsigned char tempBuffer[8];
			DWORD bytesRead = 0;

			// 读取最多8个字节
			if (ReadByte(tempBuffer, 8)) {
				// 将新读取的数据添加到缓冲区
				memcpy(buffer + bufferPos, tempBuffer, 8);
				bufferPos += 8;
			}
		}
		//unsigned char c;
		//if (ReadByte(&c, 1)) {
		//	buffer[bufferPos++] = c;
		//}
		//// 避免 buffer 溢出
		//if (bufferPos > BUFFER_SIZE - 8) {
		//	memmove(buffer, buffer + 1, --bufferPos);
		//}

		// 检查是否有有效帧
		if (findValidFrame(buffer, bufferPos, result)) {
			// 找到有效帧，调整缓冲区以移除已处理的数据
			result = convertToSingle(buffer);
			memmove(buffer, buffer + findFrameStart(buffer, bufferPos), bufferPos - findFrameStart(buffer, bufferPos));
			bufferPos -= findFrameStart(buffer, bufferPos);
			return result;
		}

		/*if (findValidFrame(buffer, bufferPos, result)) {
			result = convertToSingle(buffer);
			DWORD start = findFrameStart(buffer, bufferPos);
			memmove(buffer, buffer + start + 8, bufferPos - (start + 8));
			bufferPos -= (start + 8);
			return result;
		}*/

		// 如果缓冲区已满，移除最旧的数据
		if (bufferPos >= BUFFER_SIZE - 8) {
			memmove(buffer, buffer + 8, bufferPos - 8);
			bufferPos -= 8;
		}

		// 避免CPU占用过高
		Sleep(1);
	}
}

// 查找缓冲区中第一个有效帧的起始位置
// 如果没有找到有效帧，返回缓冲区长度（即没有有效帧）
DWORD findFrameStart(unsigned char* buffer, DWORD bufferSize) {
	for (DWORD i = 0; i <= bufferSize - 8; i++) {
		//if ((buffer[i] & 0xFF) == 0xAA && (buffer[i + 1] & 0xFF) == 0x55 && (buffer[i + 2] & 0xFF) == 0x11) {
		if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2]== 0x11) {
			//if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2] == 0x11) {
			return i;
		}
	}
	return bufferSize; // 没有找到有效帧
}

